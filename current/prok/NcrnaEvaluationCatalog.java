package prok;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
import java.util.Locale;

import map.ObjectIntMap;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.ConfigInput;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import structures.ByteBuilder;
import structures.IntList;

/** Joins evaluation metadata without constructing biological reference annotations.
 * Forked from template.A_SampleByteFile; one metadata thread owns all state.
 * Only equal resolved filesystem paths are merged. Assembly claims identify
 * possible aliases, never automatic sequence equivalence; TaxID is metadata only.
 * No genome sequence is read or hashed. Every supplied membership is preserved.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaEvaluationCatalog{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final NcrnaEvaluationCatalog tool=new NcrnaEvaluationCatalog(args);
		try{tool.process(timer);}
		finally{Shared.closeStream(tool.outstream);}
	}

	/** Parses standard BBTools options; out is a new directory, never an overwrite. */
	public NcrnaEvaluationCatalog(String[] args){
		final PreParser pp=new PreParser(args, getClass(), false);
		args=pp.args;
		outstream=pp.outstream;
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				String b=split.length>1 ? split[1] : null;
				if(b!=null && b.equalsIgnoreCase("null")){b=null;}
				if(a.equals("requirefiles")){requireFiles=Parse.parseBoolean(b);}
				else if(a.equals("resolvealiases")){resolveAliases=Parse.parseBoolean(b);}
				else if(parser.parse(arg, a, b)){
					// Standard input, output and JVM-independent BBTools options.
				}else{throw new IllegalArgumentException("Unknown parameter: "+arg);}
			}
			in=parser.in1;
			out=parser.out1;
			if(in==null || out==null || new File(out).exists()){
				throw new IllegalArgumentException("Required: in=source_specs.tsv out=new_directory");
			}
			if(parser.overwrite || parser.append || parser.in2!=null || parser.out2!=null
					|| parser.maxReads>=0 || parser.samplerate!=1){
				throw new IllegalArgumentException("A catalog requires complete unpaired metadata and fresh outputs; no overwrite, append or sampling");
			}
		}catch(RuntimeException | Error failure){
			Shared.closeStream(outstream);
			throw failure;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------          Catalog join        ----------------*/
	/*--------------------------------------------------------------*/

	/** Reads all metadata before creating outputs; missing payloads remain inspectable. */
	private void process(final Timer timer){
		assert(members.isEmpty() && genomes.isEmpty()) : "Each catalog instance owns one complete join";
		readSources();
		for(final Source source:sources){readMembers(source);}
		if(members.isEmpty()){throw new IllegalArgumentException("The supplied sources contain no memberships");}
		assert(byPath.size()==genomes.size()) : "Every resolved path must have exactly one output row";
		if(byPath.contains(normalPath(new File(out)))){
			throw new IllegalArgumentException("Output directory collides with an input FASTA path: "+out);
		}
		if(!new File(out).mkdir()){throw new IllegalArgumentException("Cannot create fresh output directory: "+out);}
		writeCatalog();
		writeAliases();
		if(resolveAliases){writeAliasResolution();}
		writeSummary();
		timer.stop();
		outstream.println("Catalog: sources="+sources.size()+" memberships="+members.size()
			+" unique_fasta_paths="+genomes.size()+" unavailable_paths="+unavailable+" time="+timer);
		if(requireFiles && unavailable>0){
			throw new IllegalStateException("Catalog contains "+unavailable+" unavailable FASTA paths; inspect genomes.tsv before scanning");
		}
	}

	/** Each specification supplies its own schema and record count, not an inferred panel size. */
	private void readSources(){
		assert(sources.isEmpty()) : "Source order establishes reproducible first-seen path identifiers";
		final File parent=new File(in).getAbsoluteFile().getParentFile();
		final ObjectSet<String> names=new ObjectSet<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(in)){
			final byte[] header=input.reader.nextLine();
			if(header==null){throw new IllegalArgumentException("Empty source specification: "+in);}
			p.set(header);
			final int width=p.terms();
			if((width!=5 && width!=6) || !p.termEquals("source", 0) || !p.termEquals("format", 1)
					|| !p.termEquals("compartment", 2) || !p.termEquals("path", 3) || !p.termEquals("expected_records", 4)
					|| (width==6 && !p.termEquals("digest_kind", 5))){
				throw new IllegalArgumentException("Expected source, format, compartment, path, expected_records TSV header, optionally followed by digest_kind");
			}
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				if(row.length==0){continue;}
				p.set(row);
				if(p.terms()!=width){throw new IllegalArgumentException("Source specification width changed");}
				final String name=p.parseString(0), format=p.parseString(1), compartment=p.parseString(2);
				final String path=p.parseString(3);
				final int expected=p.parseInt(4);
				final String digestKind=width==6 ? p.parseString(5) : "unknown";
				if(!digestKind.equals("file_bytes") && !digestKind.equals("unknown")){
					throw new IllegalArgumentException("digest_kind must be file_bytes or unknown: "+name);
				}
				if(name.isEmpty() || !names.add(name) || compartment.isEmpty() || path.isEmpty() || expected<0){
					throw new IllegalArgumentException("Source requires a unique name, compartment, path and nonnegative count: "+name);
				}
				if(!format.equals("header") && !format.equals("config") && !format.equals("taxonomy7")){
					throw new IllegalArgumentException("Unsupported source format: "+format);
				}
				sources.add(new Source(name, format, compartment, relative(parent, path), expected, digestKind));
			}
		}
		if(sources.isEmpty()){throw new IllegalArgumentException("No source specifications were supplied");}
	}

	/** Headered rosters, direct caller inputs and the documented seven-field taxonomy format. */
	private void readMembers(final Source source){
		assert(source.expected>=0) : "readSources validates the complete per-source denominator";
		final int before=members.size();
		if(source.format.equals("config")){
			try(ConfigInput input=new ConfigInput(source.file.getPath())){
				while(input.next()){
					if(input.key.equalsIgnoreCase("config")){
						throw new IllegalArgumentException("Nested configuration inputs must be listed explicitly: "+source.file);
					}
					if(!input.key.equalsIgnoreCase("in")){continue;}
					if(input.value.indexOf(',')>=0){
						throw new IllegalArgumentException("Use one FASTA per in= line, not a comma list: "+source.file);
					}
					addMember(source, input.value, NA, NA, NA, source.compartment, NA, NA);
				}
			}
		}else{
			final LineParser1 p=new LineParser1('\t');
			try(Input input=new Input(source.file.getPath())){
				final Columns columns=source.format.equals("header") ? new Columns(input.reader.nextLine(), source.file) : null;
				for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
					if(row.length==0){continue;}
					p.set(row);
					if(columns!=null){
						if(p.terms()!=columns.width){throw new IllegalArgumentException("Roster width changed: "+source.file);}
						addMember(source, p.parseString(columns.fasta), field(p, columns.genome), field(p, columns.assembly),
							field(p, columns.taxid), columns.compartment<0 ? source.compartment : field(p, columns.compartment),
							field(p, columns.gff), field(p, columns.digest));
					}else{
						if(p.terms()!=7){throw new IllegalArgumentException("taxonomy7 needs seven headerless fields: "+source.file);}
						addMember(source, p.parseString(5), field(p, 6), NA, field(p, 2), source.compartment, NA, NA);
					}
				}
			}
		}
		if(members.size()-before!=source.expected){
			throw new IllegalArgumentException("Source "+source.name+" has "+(members.size()-before)+" records; expected "+source.expected);
		}
	}

	/** Resolves file identity once per supplied absolute path; every metadata row survives. */
	private void addMember(final Source source, final String fasta, final String sourceGenome, final String assembly,
		final String taxid, final String compartment, final String gff, final String digest){
		assert(source.file.isAbsolute()) : "Relative FASTA paths are anchored to the explicit source file";
		if(fasta.isEmpty() || fasta.equals(NA)){throw new IllegalArgumentException("Missing FASTA in source "+source.name);}
		if(!digest.equals(NA) && !isSha80(digest)){
			throw new IllegalArgumentException("Reported provenance must be sha80 in source "+source.name);
		}
		final File file=relative(source.file.getParentFile(), fasta);
		// Keep the supplied spelling: link/../file cannot be simplified before
		// resolving link, because its parent may be outside the apparent directory.
		final String supplied=file.getPath();
		Genome genome=suppliedPaths.get(supplied);
		if(genome==null){
			final Genome candidate=inspect(file);
			genome=byPath.get(candidate.path);
			if(genome==null){
				genome=candidate;
				genome.id=genomes.size()+1;
				byPath.put(genome.path, genome);
				genomes.add(genome);
				if(!genome.state.equals("PRESENT")){unavailable++;}
			}
			suppliedPaths.put(supplied, genome);
		}
		genome.members++;
		if(!digest.equals(NA)){
			genome.unknownDigestKind|=!source.digestKind.equals("file_bytes");
			if(genome.digest==null){genome.digest=digest;}
			else if(!genome.digest.equals(digest)){genome.conflictingDigests=true;}
		}
		members.add(new Member(source, ++source.count, genome, fasta, sourceGenome, assembly, taxid, compartment, gff, digest));
		if(isAccession(assembly)){
			if(!genome.claims.contains(assembly)){genome.claims.add(assembly);}
			IntList ids=assemblyPaths.get(assembly);
			if(ids==null){
				ids=new IntList();
				assemblyPaths.put(assembly, ids);
				assemblyOrder.add(assembly);
			}
			if(!ids.contains(genome.id)){ids.add(genome.id);}
		}
	}

	/** Filesystem metadata only: PRESENT does not establish valid or equivalent sequence content. */
	private static Genome inspect(final File file){
		assert(file.isAbsolute()) : "A file identity cannot depend on the process working directory";
		Path path=file.toPath();
		String state;
		long bytes=-1;
		try{
			path=path.toRealPath();
			final BasicFileAttributes attrs=Files.readAttributes(path, BasicFileAttributes.class);
			bytes=attrs.size();
			state=!attrs.isRegularFile() ? "NOT_REGULAR" : !Files.isReadable(path) ? "UNREADABLE" : bytes==0 ? "EMPTY" : "PRESENT";
		}catch(NoSuchFileException failure){
			state="MISSING";
		}catch(IOException failure){
			state="UNRESOLVED_PATH";
		}
		return new Genome(path.toString(), state, bytes);
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	/** Writes distinct file paths and the complete many-to-one source membership mapping. */
	private void writeCatalog(){
		assert(!members.isEmpty()) : "A completed source join must precede catalog output";
		final ByteBuilder b=new ByteBuilder();
		try(Output paths=new Output(out+"/genomes.tsv"); Output rows=new Output(out+"/memberships.tsv")){
			paths.writer.println("input_id\tresolved_or_absolute_fasta\tpath_state\tfile_bytes\tmemberships");
			for(final Genome g:genomes){
				paths.writer.print(b.clear().append(g.id).tab().append(g.path).tab().append(g.state).tab().append(g.bytes).tab().append(g.members).nl());
			}
			rows.writer.println("source\tmember_index\tinput_id\tsupplied_fasta\tsource_genome\tsource_assembly\ttaxid\tdeclared_compartment\treported_gff\treported_sha80");
			for(final Member m:members){
				rows.writer.print(b.clear().append(m.source.name).tab().append(m.index).tab().append(m.genome.id).tab()
					.append(m.fasta).tab().append(m.sourceGenome).tab().append(m.assembly).tab().append(m.taxid).tab()
					.append(m.compartment).tab().append(m.gff).tab().append(m.digest).nl());
			}
		}
		try(Output output=new Output(out+"/sources.tsv")){
			output.writer.println("source\tformat\tcompartment\tpath\texpected_records\tobserved_records"
				+(resolveAliases ? "\tdigest_kind" : ""));
			for(final Source source:sources){
				b.clear().append(source.name).tab().append(source.format).tab().append(source.compartment)
					.tab().append(source.file.getPath()).tab().append(source.expected).tab().append(source.count);
				if(resolveAliases){b.tab().append(source.digestKind);}
				output.writer.print(b.nl());
			}
		}
	}

	/** Lists candidates for later alias adjudication, without coalescing assembly or TaxID claims. */
	private void writeAliases(){
		assert(aliasGroups==0 && conflictGroups==0) : "Alias reporting is a single pass over the completed join";
		final ByteBuilder b=new ByteBuilder();
		try(Output output=new Output(out+"/aliases.tsv")){
			output.writer.println("issue\tassembly_claim\tinput_id\tfasta");
			for(final String assembly:assemblyOrder){
				final IntList ids=assemblyPaths.get(assembly);
				if(ids.size<2){continue;}
				aliasGroups++;
				for(int i=0; i<ids.size; i++){
					final Genome g=genomes.get(ids.get(i)-1);
					output.writer.print(b.clear().append("ASSEMBLY_MULTIPLE_PATHS").tab().append(assembly).tab().append(g.id).tab().append(g.path).nl());
				}
			}
			for(final Genome g:genomes){
				if(g.claims.size()<2){continue;}
				conflictGroups++;
				for(final String assembly:g.claims){
					output.writer.print(b.clear().append("PATH_MULTIPLE_ASSEMBLIES").tab().append(assembly).tab().append(g.id).tab().append(g.path).nl());
				}
			}
		}
	}

	/** Reports one live-path candidate only when all participating paths have matching
	 * recorded identities and no conflicting claims. Does not merge rows, verify
	 * sequence contents, or relax requireFiles; unresolved groups retain a reason.
	 */
	private void writeAliasResolution(){
		assert(!members.isEmpty()) : "Alias choices require all source claims and digest disagreements to be collected";
		final ByteBuilder b=new ByteBuilder();
		try(Output output=new Output(out+"/alias_resolution.tsv")){
			output.writer.println("assembly_claim\tinput_id\toriginal_fasta\tcandidate_input_id\tcandidate_fasta\trecorded_sha80\tstatus");
			for(final String assembly:assemblyOrder){
				final IntList ids=assemblyPaths.get(assembly);
				if(ids.size<2){continue;}
				String digest=null;
				Genome candidate=null;
				int live=0;
				boolean missing=false, conflicting=false, different=false, multipleClaims=false, unknownKind=false;
				for(int i=0; i<ids.size; i++){
					final Genome g=genomes.get(ids.get(i)-1);
					multipleClaims|=g.claims.size()!=1;
					conflicting|=g.conflictingDigests;
					unknownKind|=g.unknownDigestKind;
					missing|=g.digest==null;
					if(digest==null){digest=g.digest;}
					else if(g.digest!=null && !digest.equals(g.digest)){different=true;}
					if(g.state.equals("PRESENT")){
						candidate=g;
						live++;
					}
				}
				final String status;
				if(multipleClaims){status="CONFLICTING_ASSEMBLY_CLAIMS";}
				else if(conflicting){status="CONFLICTING_PATH_DIGESTS";}
				else if(missing){status="MISSING_RECORDED_DIGEST";}
				else if(unknownKind){status="UNDECLARED_DIGEST_KIND";}
				else if(different){status="DIFFERENT_RECORDED_DIGESTS";}
				else if(live==0){status="NO_LIVE_PATH";}
				else if(live>1){status="MULTIPLE_LIVE_PATHS";}
				else{status="MATCHING_METADATA_ONE_LIVE_PATH";}
				final boolean matched=status.equals("MATCHING_METADATA_ONE_LIVE_PATH");
				assert(!matched || candidate!=null && digest!=null) : "A metadata match requires one live path and a shared recorded identity";
				for(int i=0; i<ids.size; i++){
					final Genome g=genomes.get(ids.get(i)-1);
					output.writer.print(b.clear().append(assembly).tab().append(g.id).tab().append(g.path).tab()
						.append(matched ? candidate.id : -1).tab().append(matched ? candidate.path : NA).tab()
						.append(matched ? digest : NA).tab().append(status).nl());
				}
			}
		}
	}

	/** Conserves memberships separately from physical files and unresolved biological aliases. */
	private void writeSummary(){
		long accounted=0;
		for(final Genome g:genomes){accounted+=g.members;}
		assert(accounted==members.size()) : "Every input membership must be attached to exactly one file row: "+accounted+" versus "+members.size();
		try(Output output=new Output(out+"/summary.tsv")){
			output.writer.println("sources\tmemberships\tunique_fasta_paths\tunavailable_paths\tassembly_alias_groups\tconflicting_path_groups\tsequence_equivalence");
			output.writer.print(new ByteBuilder().append(sources.size()).tab().append(members.size()).tab().append(genomes.size())
				.tab().append(unavailable).tab().append(aliasGroups).tab().append(conflictGroups).append("\tNOT_TESTED\n"));
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Helpers           ----------------*/
	/*--------------------------------------------------------------*/

	private static String field(final LineParser1 p, final int index){
		if(index<0){return NA;}
		final String value=p.parseString(index);
		return value.isEmpty() ? NA : value;
	}

	private static File relative(final File parent, final String value){
		final File file=new File(value);
		return (file.isAbsolute() ? file : new File(parent, value)).getAbsoluteFile();
	}

	private static String normalPath(final File file){
		try{return file.getCanonicalPath();}
		catch(IOException failure){throw new IllegalArgumentException("Cannot resolve output path: "+file, failure);}
	}

	/** Only explicit versioned accessions produce alias candidates; names and TaxIDs do not. */
	private static boolean isAccession(final String value){
		assert(value!=null) : "Missing metadata is represented by field() as NA before alias classification";
		int i=value.startsWith("GCF_") || value.startsWith("GCA_") ? 4 : value.startsWith("NC_") ? 3 : -1;
		if(i<0){return false;}
		final int dot=value.indexOf('.', i);
		if(dot<=i || dot==value.length()-1){return false;}
		for(; i<value.length(); i++){
			if(i!=dot && (value.charAt(i)<'0' || value.charAt(i)>'9')){return false;}
		}
		return true;
	}

	private static boolean isSha80(final String value){
		assert(value!=null) : "Digest fields are normalized before provenance validation";
		if(value.length()!=20){return false;}
		for(int i=0; i<value.length(); i++){
			final char c=value.charAt(i);
			if(!(c>='0' && c<='9' || c>='a' && c<='f')){return false;}
		}
		return true;
	}

	/** Maps retained metadata fields once per header; unknown extra columns are harmless. */
	private static final class Columns{
		Columns(final byte[] header, final File file){
			if(header==null){throw new IllegalArgumentException("Missing roster header: "+file);}
			final LineParser1 p=new LineParser1('\t');
			p.set(header);
			width=p.terms();
			final ObjectIntMap<String> fields=new ObjectIntMap<String>(String.class);
			for(int i=0; i<width; i++){
				final String name=p.parseString(i);
				if(name.isEmpty() || fields.contains(name)){throw new IllegalArgumentException("Duplicate or empty roster header in "+file);}
				fields.put(name, i);
			}
			fasta=fields.get("fasta");
			if(fasta<0){throw new IllegalArgumentException("Roster requires a named fasta column: "+file);}
			genome=fields.get("genome");
			assembly=fields.get("assembly");
			taxid=fields.get("taxid");
			compartment=fields.get("compartment");
			gff=fields.get("gff");
			int hash=-1;
			for(final String name:new String[]{"sha80", "fasta_sha80", "selected_fasta_sha80"}){
				final int index=fields.get(name);
				if(index<0){continue;}
				if(hash>=0){throw new IllegalArgumentException("Ambiguous FASTA digest columns in "+file);}
				hash=index;
			}
			digest=hash;
			assert(fasta>=0 && fasta<width) : "Roster rows use the header's FASTA index for byte-field extraction";
		}
		final int width, fasta, genome, assembly, taxid, compartment, gff, digest;
	}

	/** One explicitly named panel input, retaining its own record denominator. */
	private static final class Source{
		Source(String name_, String format_, String compartment_, File file_, int expected_, String digestKind_){
			name=name_;
			format=format_;
			compartment=compartment_;
			file=file_;
			expected=expected_;
			digestKind=digestKind_;
		}
		final String name, format, compartment, digestKind;
		final File file;
		final int expected;
		int count;
	}

	/** One resolved path, with reported assembly claims retained for later review. */
	private static final class Genome{
		Genome(String path_, String state_, long bytes_){
			path=path_;
			state=state_;
			bytes=bytes_;
		}
		final String path, state;
		final long bytes;
		final ArrayList<String> claims=new ArrayList<String>();
		/** First explicit source digest; disagreement in any later membership blocks resolution. */
		String digest;
		boolean conflictingDigests;
		boolean unknownDigestKind;
		int id, members;
	}

	/** One retained source row; metadata is reported without upgrading it to verified biology. */
	private static final class Member{
		Member(Source source_, int index_, Genome genome_, String fasta_, String sourceGenome_, String assembly_,
			String taxid_, String compartment_, String gff_, String digest_){
			source=source_;
			index=index_;
			genome=genome_;
			fasta=fasta_;
			sourceGenome=sourceGenome_;
			assembly=assembly_;
			taxid=taxid_;
			compartment=compartment_;
			gff=gff_;
			digest=digest_;
		}
		final Source source;
		final Genome genome;
		final int index;
		final String fasta, sourceGenome, assembly, taxid, compartment, gff, digest;
	}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	private String in, out;
	private PrintStream outstream=System.err;
	private boolean requireFiles=true;
	/** Reports metadata-supported choices only; never changes membership or availability guards. */
	private boolean resolveAliases=false;
	private int unavailable, aliasGroups, conflictGroups;
	private final ArrayList<Source> sources=new ArrayList<Source>();
	private final ArrayList<Genome> genomes=new ArrayList<Genome>();
	private final ArrayList<Member> members=new ArrayList<Member>();
	private final ObjectMap<String,Genome> byPath=new ObjectMap<String,Genome>(String.class, Genome.class);
	private final ObjectMap<String,Genome> suppliedPaths=new ObjectMap<String,Genome>(String.class, Genome.class);
	private final ObjectMap<String,IntList> assemblyPaths=new ObjectMap<String,IntList>(String.class, IntList.class);
	private final ArrayList<String> assemblyOrder=new ArrayList<String>();
	private static final String NA="NA";
}
