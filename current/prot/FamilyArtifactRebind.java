package prot;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.LinkOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

import ddl.DDLFamilySearch.StructuralInputs;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;

import parse.LineParser1;
import parse.Parse;
import parse.PreParser;

import prok.StopTailCascadeResources;

import prot.FamilyCalibration.VerificationMode;

import shared.Shared;

import structures.ByteBuilder;

/** Explicit byte-preserving path edits in legacy classifier/cascade text artifacts.
 * Family data and fitted values are never rewritten or loaded.
 * @author Sayu
 */
public final class FamilyArtifactRebind{
	/** Parse only explicit edit-table/output bindings; close the diagnostic stream on every exit. */
	public static void main(String[] args)throws Exception{
		final PreParser pp=new PreParser(args, FamilyArtifactRebind.class, false);
		try{
			Path in=null;
			Path out=null;
			for(String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				require(pair.length==2 && !pair[1].isEmpty(), "Require in=bindings.tsv out=fresh_directory");
				if(pair[0].equalsIgnoreCase("in") && in==null){
					in=Paths.get(pair[1]);
				}else if(pair[0].equalsIgnoreCase("out") && out==null){
					out=Paths.get(pair[1]);
				}else{
					throw new IllegalArgumentException("Unknown or duplicate option: "+pair[0]);
				}
			}
			require(in!=null && out!=null, "Require in=bindings.tsv out=fresh_directory");
			Shared.setThreads(1);
			rebind(in, out);
		}finally{
			Shared.closeStream(pp.outstream);
		}
	}

	/** All source graph checks precede candidate creation; PASS certifies only a complete batch. */
	public static void rebind(Path bindings, Path output)throws Exception{
		require(output.isAbsolute() && !Files.exists(output, LinkOption.NOFOLLOW_LINKS), "Fresh absolute output directory required");
		final Path out=output.normalize().getParent().toRealPath().resolve(output.getFileName());
		require(!Files.exists(out, LinkOption.NOFOLLOW_LINKS), "Output already exists");
		final StructuralInputs stable=new StructuralInputs();
		stable.track(regular(bindings));
		final LinkedHashMap<String, Artifact> artifacts=bindings(bindings, out);
		final HashMap<Path, Artifact> bySource=new HashMap<Path, Artifact>();
		final HashSet<Path> targets=new HashSet<Path>();
		final HashMap<Path, Document> documents=new HashMap<Path, Document>();
		for(Artifact a:artifacts.values()){
			require(bySource.put(a.source, a)==null && targets.add(a.target), "Duplicate artifact source or target");
			a.before=document(a.source, a.cascade, stable);
			documents.put(a.source, a.before);
			for(Map.Entry<String, Edit> entry:a.edits.entrySet()){
				final Field field=a.before.fields.get(entry.getKey());
				require(field!=null && field.value.equals(entry.getValue().expected), "Missing key or wrong old value: "+a.id+"/"+entry.getKey());
			}
		}
		// Inspect the complete original text graph BEFORE invoking any structural reader.
		final HashSet<Path> referenced=new HashSet<Path>();
		for(Artifact a:artifacts.values()){
			preflight(a.before, documents, stable, new HashSet<Path>(), referenced);
		}
		final ArrayList<Artifact> roots=new ArrayList<Artifact>();
		for(Artifact a:artifacts.values()){
			if(!a.cascade){
				a.classifier=new FamilyScoringPackage(a.source, VerificationMode.STRUCTURAL);
			}else if(!referenced.contains(a.source)){
				a.resources=StopTailCascadeResources.readStructural(a.source.toString(), "unmeasured");
				a.roles=a.resources.resolvedResourceBindings();
				roots.add(a);
			}
		}
		final HashSet<Path> covered=new HashSet<Path>();
		for(Artifact root:roots){
			for(Map.Entry<String, Path> role:root.roles.entrySet()){
				if(manifest(role.getKey())){
					covered.add(role.getValue().toRealPath());
				}
				protect(out, role.getValue());
			}
		}
		for(Artifact a:artifacts.values()){
			if(a.cascade){
				require(covered.contains(a.source), "Cascade edit needs its composite root in this batch: "+a.id);
			}else{
				protect(out, a.source);
				for(String value:a.classifier.resolvedOptions().values()){
					if(Paths.get(value).isAbsolute() && Files.exists(Paths.get(value))){
						protect(out, Paths.get(value));
					}
				}
			}
		}
		stable.verify();
		Files.createDirectory(out);
		for(Artifact a:artifacts.values()){
			write(a.target, a.before.replace(a.edits));
		}
		// Candidate children are now present at their declared final paths. Repeat the
		// no-hash preflight before binding that graph, including unchanged external children.
		final HashMap<Path, Document> candidates=new HashMap<Path, Document>();
		for(Artifact a:artifacts.values()){
			a.after=document(a.target, a.cascade, stable);
			candidates.put(a.target, a.after);
			preflight(a.after, candidates, stable, new HashSet<Path>(), new HashSet<Path>());
		}
		for(Artifact a:artifacts.values()){
			if(!a.cascade){
				final FamilyScoringPackage next=new FamilyScoringPackage(a.target, VerificationMode.STRUCTURAL);
				final HashMap<String, String> expected=a.classifier.resolvedOptions();
				final Edit edit=a.edits.get("blacklist");
				if(edit!=null){
					expected.put("blacklist", Paths.get(edit.replacement).normalize().toString());
				}
				require(expected.equals(next.resolvedOptions()), "Unedited classifier binding changed after relocation: "+a.id);
				a.classifier.verify();
				next.verify();
			}
		}
		final HashSet<Path> candidateCovered=new HashSet<Path>();
		for(Artifact root:roots){
			final StopTailCascadeResources next=StopTailCascadeResources.readStructural(root.target.toString(), "unmeasured");
			final Map<String, Path> actual=next.resolvedResourceBindings();
			for(Map.Entry<String, Path> entry:actual.entrySet()){
				if(manifest(entry.getKey())){
					candidateCovered.add(entry.getValue().toRealPath());
				}
			}
			require(root.roles.keySet().equals(actual.keySet()), "Resolved cascade role roster changed");
			for(Map.Entry<String, Path> role:root.roles.entrySet()){
				Path expected=role.getValue();
				final String key=role.getKey();
				if(manifest(key)){
					final Artifact moved=bySource.get(expected.toRealPath());
					if(moved!=null){
						expected=moved.target;
					}
					if(key.endsWith("short40.manifest")){
						final String parent=key.substring(0, key.length()-"short40.manifest".length())+"manifest";
						final Artifact owner=bySource.get(root.roles.get(parent).toRealPath());
						if(owner!=null && owner.edits.containsKey("short40manifest")){
							expected=Paths.get(owner.edits.get("short40manifest").replacement);
						}
					}
				}else if(key.equals("blacklist") || key.endsWith(".blacklist")){
					final String owner=key.substring(0, key.length()-"blacklist".length())+"manifest";
					final Artifact moved=bySource.get(root.roles.get(owner).toRealPath());
					if(moved!=null && moved.edits.containsKey("blacklist")){
						expected=Paths.get(moved.edits.get("blacklist").replacement);
					}
				}
				require(Files.isSameFile(expected, actual.get(key)), "Unedited cascade resource changed after relocation: "+key);
				if(manifest(key)){
					final Document original=documents.get(role.getValue().toRealPath());
					final Artifact declared=bySource.get(role.getValue().toRealPath());
					final byte[] wanted=declared==null?original.bytes:original.replace(declared.edits).toBytes();
					require(java.util.Arrays.equals(wanted, Files.readAllBytes(actual.get(key))), "Undeclared cascade text/policy change: "+key);
				}
			}
			root.resources.verify();
			next.verify();
		}
		for(Artifact a:artifacts.values()){
			if(a.cascade){
				require(candidateCovered.contains(a.target.toRealPath()), "Rebound cascade is not covered by its candidate composite graph: "+a.id);
			}
		}
		stable.verify();
		final ByteBuilder receipt=new ByteBuilder(HEADER).nl().append(COLUMNS).nl();
		int count=0;
		for(Artifact a:artifacts.values()){
			for(Map.Entry<String, Edit> entry:a.edits.entrySet()){
				final Edit edit=entry.getValue();
				count++;
				receipt.append(a.id).tab().append(a.cascade?"cascade":"classifier").tab().append(utf8(a.source.toString())).tab()
				.append(utf8(a.target.getFileName().toString())).tab().append(entry.getKey()).tab().append(utf8(edit.expected)).tab().append(utf8(edit.replacement)).nl();
			}
		}
		write(out.resolve("bindings.tsv"), receipt);
		stable.verify();
		write(out.resolve("PASS"), new ByteBuilder("FAMILY_ARTIFACT_REBIND_PASS artifacts=").append(artifacts.size())
		.append(" changed_path_fields=").append(count).append(" classifier_runs=0\n"));
	}

	/** Reuse a borrowed row parser; retain only validated edit descriptors grouped by artifact ID. */
	private static LinkedHashMap<String, Artifact> bindings(Path input, Path out)throws Exception{
		final LinkedHashMap<String, Artifact> artifacts=new LinkedHashMap<String, Artifact>();
		final ByteFile in=ByteFile.makeByteFile(input.toString(), false);
		final LineParser1 parser=new LineParser1('\t');
		try{
			require(java.util.Arrays.equals(utf8(HEADER), in.nextLine()) && java.util.Arrays.equals(utf8(COLUMNS), in.nextLine()), "Expected family rebind format/version and columns");
			for(byte[] line; (line=in.nextLine())!=null; ){
				parser.set(line);
				require(parser.terms()==7, "Malformed rebind row");
				final String id=text(parser, 0);
				final String name=text(parser, 3);
				require(safe(id) && safe(name) && !name.equals("PASS") && !name.equals("bindings.tsv"), "Unsafe artifact/target");
				final boolean cascade=parser.termEquals("cascade", 1);
				require(cascade || parser.termEquals("classifier", 1), "Unknown artifact format");
				require(parser.termEquals("blacklist", 4) || cascade && parser.termEquals("short40manifest", 4), "Only declared path keys may be rebound");
				final String key=parser.termEquals("blacklist", 4)?"blacklist":"short40manifest";
				final String expected=text(parser, 5);
				final String replacement=text(parser, 6);
				final Path declared=Paths.get(text(parser, 2));
				require(!expected.isEmpty() && declared.isAbsolute() && Paths.get(replacement).isAbsolute(), "Require absolute source/replacement and exact old value");
				final Path source=regular(declared).toRealPath();
				final Path target=out.resolve(name);
				Artifact a=artifacts.get(id);
				if(a==null){
					a=new Artifact(id, source, target, cascade);
					artifacts.put(a.id, a);
				}
				require(a.source.equals(source) && a.target.equals(target) && a.cascade==cascade, "Conflicting artifact rows");
				require(a.edits.put(key, new Edit(expected, replacement))==null, "Duplicate requested key");
			}
		}finally{
			require(!in.close(), "Failed reading rebind table");
		}
		require(!artifacts.isEmpty(), "Empty rebind batch");
		return artifacts;
	}

	/** Parse only byte spans and child references; semantic schemas belong to existing readers. */
	private static Document document(Path file, boolean cascade, StructuralInputs stable)throws Exception{
		regular(file);
		stable.track(file);
		final byte[] bytes=Files.readAllBytes(file);
		final LinkedHashMap<String, Field> fields=new LinkedHashMap<String, Field>();
		int start=0;
		int number=0;
		while(start<bytes.length){
			int end=start;
			while(end<bytes.length && bytes[end]!='\n'){
				end++;
			}
			int lo=start;
			int hi=end;
			if(cascade){
				// ByteFile removes exactly one CR before LF; an unterminated final CR is data.
				if(end<bytes.length && hi>lo && bytes[hi-1]=='\r'){
					hi--;
				}
			}else{
				while(lo<hi && (bytes[lo]&255)<=32){
					lo++;
				}
				while(hi>lo && (bytes[hi-1]&255)<=32){
					hi--;
				}
			}
			final String line=new String(bytes, lo, hi-lo, StandardCharsets.UTF_8);
			if(number++==0 && cascade){
				require(line.matches("#bbgap_tail_cascade\\t[1-8]"), "Unsupported cascade header (fitted forms are not rebound)");
			}
			if(lo<hi && bytes[lo]!='#'){
				int sep=lo;
				while(sep<hi && bytes[sep]!=(cascade?'\t':'=')){
					sep++;
				}
				require(sep>lo && sep<hi, "Malformed artifact field: "+file);
				final String rawKey=new String(bytes, lo, sep-lo, StandardCharsets.UTF_8);
				final String key=(cascade?rawKey:rawKey.trim()).toLowerCase(Locale.ROOT);
				int valueStart=sep+1;
				if(!cascade){
					while(valueStart<hi && (bytes[valueStart]&255)<=32){
						valueStart++;
					}
				}
				final String value=new String(bytes, valueStart, hi-valueStart, StandardCharsets.UTF_8);
				require(!value.isEmpty() && !value.contains("\t") && fields.put(key, new Field(value, valueStart, hi))==null, "Empty/duplicate artifact key or extra column");
				require(!key.equals("uncalibratedmatrixfile") && !key.equals("uncalibratedmatrixfilesha80") &&
				!key.equals("familymeans") && !key.equals("familymeanssha80"), "Unsupported diagnostic-matrix or scalar-mean binding before structural reader");
			}
			start=end+1;
		}
		require(number>0, "Empty artifact");
		stable.verify();
		return new Document(file, cascade, bytes, fields);
	}

	/** Walk text children before reader entry, tracking stability and rejecting recursion or unsafe forms. */
	private static void preflight(Document doc, Map<Path, Document> documents, StructuralInputs stable,
	Set<Path> visiting, Set<Path> referenced)throws Exception{
		require(visiting.add(doc.path), "Cyclic cascade manifest graph");
		if(doc.cascade){
			for(String key:new String[]{
					"short40manifest", "short100manifest", "short150manifest"
			}){
				final Field field=doc.fields.get(key);
				if(field==null){
					continue;
				}
				final Path child=regular(doc.path.getParent().resolve(field.value).normalize()).toRealPath();
				referenced.add(child);
				Document nested=documents.get(child);
				if(nested==null){
					nested=document(child, true, stable);
					documents.put(child, nested);
				}
				require(nested.cascade, "Cascade child must be a cascade artifact");
				preflight(nested, documents, stable, visiting, referenced);
			}
		}
		visiting.remove(doc.path);
	}

	private static boolean manifest(String role){
		return role.equals("manifest") || role.endsWith(".manifest");
	}

	private static boolean safe(String value){
		return value.matches("[A-Za-z0-9][A-Za-z0-9_.-]*");
	}

	/** Decode a selected descriptor field as UTF-8; the parser itself remains borrowed. */
	private static String text(LineParser1 parser, int term){
		parser.setBounds(term);
		return new String(parser.line(), parser.a(), parser.b()-parser.a(), StandardCharsets.UTF_8);
	}

	private static byte[] utf8(String text){
		return text.getBytes(StandardCharsets.UTF_8);
	}

	/** Reject directories and leaf symlinks before any text read; return an absolute normalized path. */
	private static Path regular(Path path){
		require(Files.isRegularFile(path, LinkOption.NOFOLLOW_LINKS) && !Files.isSymbolicLink(path), "Regular nonsymlink text input required: "+path);
		return path.toAbsolutePath().normalize();
	}

	/** A fresh output root may not contain an existing resolved input. */
	private static void protect(Path out, Path input)throws Exception{
		require(!input.toRealPath().startsWith(out), "Output would contain an existing resource: "+input);
	}

	/** Wait for each fresh byte output to close successfully before publishing the batch marker. */
	private static void write(Path path, ByteBuilder bytes){
		final ByteStreamWriter out=new ByteStreamWriter(path.toString(), false, false, false);
		out.start();
		out.print(bytes);
		require(!out.poisonAndWait(), "Artifact write failed: "+path);
	}

	private static void require(boolean valid, String why){
		if(!valid){
			throw new IllegalArgumentException(why);
		}
	}

	/** Own one source/target edit group and the metadata readers that must remain stable. */
	private static final class Artifact{
		Artifact(String id_, Path source_, Path target_, boolean cascade_){
			id=id_;
			source=source_;
			target=target_;
			cascade=cascade_;
		}
		final String id;
		final Path source, target;
		final boolean cascade;
		final LinkedHashMap<String, Edit> edits=new LinkedHashMap<String, Edit>();
		Document before, after;
		FamilyScoringPackage classifier;
		StopTailCascadeResources resources;
		Map<String, Path> roles;
	}

	/** Exact old logical value and explicit replacement; never inferred from a study name. */
	private static final class Edit{
		Edit(String expected_, String replacement_){
			expected=expected_;
			replacement=replacement_;
		}
		final String expected, replacement;
	}

	/** Logical value plus its half-open byte span in the owning Document. */
	private static final class Field{
		Field(String value_, int start_, int end_){
			value=value_;
			start=start_;
			end=end_;
		}
		final String value;
		final int start, end;
	}

	/** Own immutable source bytes and physical field order; replacement copies every unedited span. */
	private static final class Document{
		Document(Path path_, boolean cascade_, byte[] bytes_, LinkedHashMap<String, Field> fields_){
			path=path_;
			cascade=cascade_;
			bytes=bytes_;
			fields=fields_;
		}
		ByteBuilder replace(Map<String, Edit> edits){
			final ByteBuilder result=new ByteBuilder(bytes.length+256);
			int from=0;
			for(Map.Entry<String, Field> entry:fields.entrySet()){
				final Edit edit=edits.get(entry.getKey());
				if(edit==null){
					continue;
				}
				final Field field=entry.getValue();
				result.append(bytes, from, field.start-from).append(edit.replacement.getBytes(StandardCharsets.UTF_8));
				from=field.end;
			}
			return result.append(bytes, from, bytes.length-from);
		}
		final Path path;
		final boolean cascade;
		final byte[] bytes;
		final LinkedHashMap<String, Field> fields;
	}

	public static final String HEADER="#bbgap_family_rebind\t1";
	public static final String COLUMNS="artifact\tformat\tsource\ttarget\tkey\texpected\treplacement";
}
