package prot;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Producer and fail-closed validator for the selection-specific reference-CDS
 * survival sidecar association.
 *
 * <p>The index key is {@code (bin_manifest_sha256, bin_id)}.  It records the
 * sidecar and its source-manifest provenance, and copies the exact native and
 * foreign source selections from the referenced bin.  The role-separated span
 * list is derived from the hash-gated source manifest's selected FASTA headers;
 * it is never accepted from the caller as proof.</p>
 */
public final class ReferenceCdsSidecarIndex {

	/** Legacy index: source span/header identifiers are stored as v1 wrote them. */
	public static final String SCHEMA_V1 = "reference_cds_sidecar_index_v1";
	/** v2 index: source span/header identifiers use IdentifierCodec percent escaping. */
	public static final String SCHEMA = "reference_cds_sidecar_index_v2";
	private static final String COLUMNS = "bin_manifest_sha256\tbin_id\tsidecar_path\tsidecar_sha256\t"
		+ "source_manifest_path\tsource_manifest_sha256\tselection_binding_sha256\t"
		+ "native_source_ids\tforeign_source_ids\tselected_source_span_ids";
	private static final String BINDING_COLUMNS = "bin_id\tsidecar_path\tsource_manifest_path";
	private static final String SOURCE_MANIFEST_COLUMNS = "role\tassembly_id\treference_gff\tassembly_fasta\tselected_shreds";
	private static final Pattern INSTANCE_SUFFIX_V1=Pattern.compile("^(.*)__instance_([A-Za-z0-9._-]+)$");
	/** v2 permits the escaped identifier bytes that can occur in an instance ID. */
	private static final Pattern INSTANCE_SUFFIX=Pattern.compile("^(.*)__instance_(.+)$");
	private static final String HEX="[0-9a-fA-F]{64}";

	private ReferenceCdsSidecarIndex() { }

	public static final class Entry {
		public final String binManifestSha256, binId, sidecarPath, sidecarSha256;
		public final String sourceManifestPath, sourceManifestSha256, selectionBindingSha256;
		public final String nativeSourceIds, foreignSourceIds, selectedSourceSpanIds;
		Entry(String binManifestSha256_, String binId_, String sidecarPath_, String sidecarSha256_,
				String sourceManifestPath_, String sourceManifestSha256_, String selectionBindingSha256_,
				String nativeSourceIds_, String foreignSourceIds_, String selectedSourceSpanIds_) {
			binManifestSha256=binManifestSha256_; binId=binId_; sidecarPath=sidecarPath_; sidecarSha256=sidecarSha256_;
			sourceManifestPath=sourceManifestPath_; sourceManifestSha256=sourceManifestSha256_;
			selectionBindingSha256=selectionBindingSha256_; nativeSourceIds=nativeSourceIds_;
			foreignSourceIds=foreignSourceIds_; selectedSourceSpanIds=selectedSourceSpanIds_;
		}
	}

	public static final class Index {
		public final String binManifestPath, binManifestSha256;
		public final int schemaVersion;
		public final boolean encoded;
		public final ArrayList<Entry> entries;
		Index(String path, String hash, int schemaVersion_, ArrayList<Entry> rows) {
			binManifestPath=path; binManifestSha256=hash; schemaVersion=schemaVersion_;
			encoded=schemaVersion_==MagQCBinManifest.SCHEMA_VERSION_ENCODED; entries=rows;
		}
	}

	private static final class Binding {
		final String binId, sidecarPath, sourceManifestPath;
		Binding(String b, String s, String m) {
			binId=b; sidecarPath=s; sourceManifestPath=m;
		}
	}

	private static final class SourceRow {
		final String role, assemblyId, referenceGff, assemblyFasta, selectedShreds;
		SourceRow(String r, String a, String g, String f, String s) {
			role=r; assemblyId=a; referenceGff=g; assemblyFasta=f; selectedShreds=s;
		}
		String key() { return role+"\t"+assemblyId; }
	}

	public static void main(String[] args) {
		if(args.length<1) { throw new RuntimeException("Usage: build binmanifest=<file> bindings=<file> out=<file> | validate index=<file>"); }
		if("build".equalsIgnoreCase(args[0])) {
			final Map<String,String> p=parseArgs(args, 1);
			build(require(p,"binmanifest"), require(p,"bindings"), require(p,"out"));
		} else if("validate".equalsIgnoreCase(args[0])) {
			final Map<String,String> p=parseArgs(args, 1);
			final Index i=load(require(p,"index"));
			System.err.println("ReferenceCdsSidecarIndex: validated " + i.entries.size()
					+ " entries; bin_manifest_sha256=" + i.binManifestSha256);
		} else {
			throw new RuntimeException("Unknown command: " + args[0]);
		}
	}

	/** Builds and immediately validates an index from a bin manifest and binding table. */
	public static String build(String binManifest, String bindingsFile, String outFile) {
		final MagQCBinManifest.Manifest manifest=MagQCBinManifest.load(binManifest);
		final String binHash=MagQCBinManifest.sha256File(binManifest);
		final HashMap<String,MagQCBinManifest.Bin> bins=new HashMap<String,MagQCBinManifest.Bin>();
		for(MagQCBinManifest.Bin b : manifest.bins) { bins.put(b.id,b); }
		final ArrayList<Binding> bindings=readBindings(bindingsFile);
		final HashSet<String> ids=new HashSet<String>();
		final ArrayList<Entry> rows=new ArrayList<Entry>();
		for(Binding b : bindings) {
			if(!ids.add(b.binId)) { throw new RuntimeException("Duplicate binding bin_id: " + b.binId); }
			final MagQCBinManifest.Bin bin=bins.get(b.binId);
			if(bin==null) { throw new RuntimeException("Binding bin_id absent from manifest: " + b.binId); }
			if(bin.schemaVersion!=manifest.schemaVersion) {
				throw new RuntimeException("Bin schema differs from manifest for " + b.binId);
			}
			final String sidecarHash=hashExisting(b.sidecarPath, "sidecar");
			final String sourceHash=hashExisting(b.sourceManifestPath, "source manifest");
			final String sidecarSourceHash=sidecarSourceManifestHash(b.sidecarPath);
			if(!sourceHash.equals(sidecarSourceHash)) {
				throw new RuntimeException("Sidecar/source-manifest hash mismatch for " + b.binId
						+ ": sidecar=" + sidecarSourceHash + " source=" + sourceHash);
			}
			final String nativeIds=bin.nativeContigs, foreignIds=bin.foreignContigs;
			final String bindingHash=selectionBindingSha256(nativeIds, foreignIds);
			final String selectedSourceSpanIds=deriveSelectedSourceSpanIds(b.sourceManifestPath, b.sidecarPath, bin);
			rows.add(new Entry(binHash,b.binId,b.sidecarPath,sidecarHash,b.sourceManifestPath,sourceHash,
					bindingHash,nativeIds,foreignIds,selectedSourceSpanIds));
		}
		write(outFile, binManifest, binHash, manifest.schemaVersion, rows);
		load(outFile);
		return MagQCBinManifest.sha256File(outFile);
	}

	/** Loads and validates every provenance, key, selection, and hash guard. */
	public static Index load(String indexFile) {
		String schema=null, binManifestPath=null, binManifestHash=null, columns=null;
		final ArrayList<Entry> rows=new ArrayList<Entry>();
		final HashSet<String> keys=new HashSet<String>();
		try(BufferedReader br=new BufferedReader(new FileReader(indexFile))) {
			for(String line; (line=br.readLine())!=null;) {
				if(line.length()==0) { continue; }
				if(line.charAt(0)=='#') {
					final String[] p=line.split("\\t",-1);
					if(line.startsWith("#schema_version\t")) {
						if(p.length!=2 || schema!=null || !(SCHEMA_V1.equals(p[1]) || SCHEMA.equals(p[1]))) { fail(indexFile,"invalid schema header"); }
						schema=p[1];
					} else if(line.startsWith("#bin_manifest\t")) {
						if(p.length!=5 || !"#".equals(p[2]) || !"sha256".equals(p[3]) || binManifestPath!=null) { fail(indexFile,"invalid bin manifest header"); }
						binManifestPath=p[1]; binManifestHash=p[4];
					} else if(line.startsWith("#columns\t")) {
						final String actualColumns=line.substring("#columns\t".length());
						if(columns!=null || !COLUMNS.equals(actualColumns)) { fail(indexFile,"invalid columns header"); }
						columns=actualColumns;
					} else {
						fail(indexFile,"unknown metadata header: "+line);
					}
					continue;
				}
				if(schema==null || binManifestPath==null || columns==null) { fail(indexFile,"data precedes required headers"); }
				final String[] p=line.split("\\t",-1);
				if(p.length!=10) { fail(indexFile,"expected 10 data fields, got "+p.length); }
				for(String s : p) { if(s.indexOf('\r')>=0 || s.indexOf('\n')>=0) { fail(indexFile,"newline in field"); } }
				final Entry e=new Entry(p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9]);
				final String key=e.binManifestSha256+"\t"+e.binId;
				if(!keys.add(key)) { fail(indexFile,"duplicate key "+key); }
				rows.add(e);
			}
		} catch(IOException ex) { throw new RuntimeException("Could not read index "+indexFile,ex); }
		if(schema==null || binManifestPath==null || binManifestHash==null || columns==null) { fail(indexFile,"incomplete headers"); }
		final String actualBinHash=hashExisting(binManifestPath,"bin manifest");
		if(!actualBinHash.equals(binManifestHash)) { fail(indexFile,"bin manifest hash mismatch"); }
		final MagQCBinManifest.Manifest manifest=MagQCBinManifest.load(binManifestPath);
		final int indexSchema=SCHEMA.equals(schema) ? MagQCBinManifest.SCHEMA_VERSION_ENCODED
			: MagQCBinManifest.SCHEMA_VERSION;
		final String expectedSchema=indexSchema==MagQCBinManifest.SCHEMA_VERSION_ENCODED ? SCHEMA : SCHEMA_V1;
		if(manifest.schemaVersion!=indexSchema || !expectedSchema.equals(schema)) {
			fail(indexFile,"index schema does not match bound bin manifest schema");
		}
		final HashMap<String,MagQCBinManifest.Bin> bins=new HashMap<String,MagQCBinManifest.Bin>();
		for(MagQCBinManifest.Bin b : manifest.bins) { bins.put(b.id,b); }
		for(Entry e : rows) {
			if(!e.binManifestSha256.equals(binManifestHash)) { fail(indexFile,"row key uses a different bin manifest hash: "+e.binId); }
			final MagQCBinManifest.Bin bin=bins.get(e.binId);
			if(bin==null) { fail(indexFile,"row bin_id absent from manifest: "+e.binId); }
			if(!bin.nativeContigs.equals(e.nativeSourceIds) || !bin.foreignContigs.equals(e.foreignSourceIds)) {
				fail(indexFile,"selection binding differs from recorded bin: "+e.binId);
			}
			final String expectedBinding=selectionBindingSha256(e.nativeSourceIds,e.foreignSourceIds);
			if(!expectedBinding.equals(e.selectionBindingSha256)) { fail(indexFile,"selection binding hash mismatch: "+e.binId); }
			final String sidecarHash=hashExisting(e.sidecarPath,"sidecar");
			if(!sidecarHash.equals(e.sidecarSha256)) { fail(indexFile,"sidecar hash mismatch: "+e.binId); }
			final String sourceHash=hashExisting(e.sourceManifestPath,"source manifest");
			if(!sourceHash.equals(e.sourceManifestSha256)) { fail(indexFile,"source manifest hash mismatch: "+e.binId); }
			if(!sourceHash.equals(sidecarSourceManifestHash(e.sidecarPath))) { fail(indexFile,"sidecar source-manifest binding mismatch: "+e.binId); }
			final String expectedSpans=deriveSelectedSourceSpanIds(e.sourceManifestPath, e.sidecarPath, bin);
			if(!expectedSpans.equals(e.selectedSourceSpanIds)) { fail(indexFile,"selected source span/header set mismatch: "+e.binId); }
		}
		return new Index(binManifestPath,binManifestHash,indexSchema,rows);
	}

	private static void write(String file, String binManifestPath, String binHash, int schemaVersion,
			ArrayList<Entry> rows) {
		try(BufferedWriter bw=new BufferedWriter(new FileWriter(file))) {
			final String schema=schemaVersion==MagQCBinManifest.SCHEMA_VERSION ? SCHEMA_V1
				: schemaVersion==MagQCBinManifest.SCHEMA_VERSION_ENCODED ? SCHEMA : null;
			if(schema==null) { throw new RuntimeException("Unrecognized sidecar index schema version: "+schemaVersion); }
			bw.write("#schema_version\t"+schema+"\n");
			bw.write("#bin_manifest\t"+binManifestPath+"\t#\tsha256\t"+binHash+"\n");
			bw.write("#columns\t"+COLUMNS+"\n");
			for(Entry e : rows) {
				bw.write(e.binManifestSha256+"\t"+e.binId+"\t"+e.sidecarPath+"\t"+e.sidecarSha256+"\t"
						+e.sourceManifestPath+"\t"+e.sourceManifestSha256+"\t"+e.selectionBindingSha256+"\t"
						+e.nativeSourceIds+"\t"+e.foreignSourceIds+"\t"+e.selectedSourceSpanIds+"\n");
			}
		} catch(IOException ex) { throw new RuntimeException("Could not write index "+file,ex); }
		try {
			final String hash=MagQCBinManifest.sha256File(file);
			Files.write(Paths.get(file+".sha256"),(hash+"  "+file+"\n").getBytes(StandardCharsets.UTF_8));
		} catch(IOException ex) { throw new RuntimeException("Could not write index hash sidecar "+file,ex); }
	}

	private static ArrayList<Binding> readBindings(String file) {
		final ArrayList<Binding> out=new ArrayList<Binding>();
		try(BufferedReader br=new BufferedReader(new FileReader(file))) {
			boolean header=false;
			for(String line; (line=br.readLine())!=null;) {
				if(line.length()==0 || line.charAt(0)=='#') { continue; }
				final String[] p=line.split("\\t",-1);
				if(!header) {
					if(p.length!=3 || !BINDING_COLUMNS.equals(line)) { throw new RuntimeException("Invalid binding columns in "+file); }
					header=true; continue;
				}
				if(p.length!=3) { throw new RuntimeException("Binding row has "+p.length+" fields in "+file); }
				for(int i=0;i<p.length;i++) { if(p[i].length()==0 || p[i].indexOf('\t')>=0) { throw new RuntimeException("Empty/malformed binding field in "+file); } }
				out.add(new Binding(p[0],p[1],p[2]));
			}
			if(!header) { throw new RuntimeException("Missing binding columns in "+file); }
		} catch(IOException ex) { throw new RuntimeException("Could not read bindings "+file,ex); }
		return out;
	}

	private static String deriveSelectedSourceSpanIds(String sourceManifestPath, String sidecarPath,
			MagQCBinManifest.Bin bin) {
		final boolean encoded=bin.schemaVersion==MagQCBinManifest.SCHEMA_VERSION_ENCODED;
		final ArrayList<SourceRow> sources=readSourceManifest(sourceManifestPath);
		validateSourceInputHashes(sourceManifestPath, sidecarPath, sources);
		final TreeSet<String> observed=new TreeSet<String>();
		for(SourceRow source : sources) {
			final File selected=resolveExistingPath(sourceManifestPath, source.selectedShreds, "selected shreds");
			for(String header : readSelectedHeaders(selected)) {
				final Matcher instance=(encoded ? INSTANCE_SUFFIX : INSTANCE_SUFFIX_V1).matcher(header);
				if(!instance.matches() || instance.group(1).length()==0) {
					throw new RuntimeException("Selected FASTA header lacks an explicit __instance_<id>: "
							+header+" in "+selected);
				}
				requireSafeSpanHeader(header, encoded);
				final String key=source.role+"|"+(encoded ? IdentifierCodec.encode(header) : header);
				if(!observed.add(key)) {
					throw new RuntimeException("Duplicate role-separated selected span/header: "+key);
				}
			}
		}
		final TreeSet<String> expected=new TreeSet<String>();
		addExpectedSelections(expected,"native",bin,false);
		addExpectedSelections(expected,"foreign",bin,true);
		if(!expected.equals(observed)) {
			throw new RuntimeException("Selected source span/header set differs for "+bin.id
					+": expected="+join(expected)+" observed="+join(observed));
		}
		return join(observed);
	}

	private static void addExpectedSelections(TreeSet<String> out, String role,
			MagQCBinManifest.Bin bin, boolean foreign) {
		final boolean encoded=bin.schemaVersion==MagQCBinManifest.SCHEMA_VERSION_ENCODED;
		for(MagQCBinManifest.Selection selection : MagQCBinManifest.decodeSelections(bin, foreign)) {
			final String header=selection.contigId+"__instance_"+selection.instanceId;
			requireSafeSpanHeader(header, encoded);
			final String key=role+"|"+(encoded ? IdentifierCodec.encode(header) : header);
			if(!out.add(key)) {
				throw new RuntimeException("Duplicate expected role-separated selection: "+key);
			}
		}
	}

	private static ArrayList<SourceRow> readSourceManifest(String file) {
		final ArrayList<SourceRow> out=new ArrayList<SourceRow>();
		final HashSet<String> keys=new HashSet<String>();
		try(BufferedReader br=new BufferedReader(new FileReader(file))) {
			boolean header=false;
			for(String line; (line=br.readLine())!=null;) {
				if(line.length()==0 || line.charAt(0)=='#') { continue; }
				final String[] p=line.split("\\t",-1);
				if(!header) {
					if(p.length!=5 || !SOURCE_MANIFEST_COLUMNS.equals(line)) {
						throw new RuntimeException("Invalid source-manifest columns in "+file);
					}
					header=true; continue;
				}
				if(p.length!=5) { throw new RuntimeException("Source-manifest row has "+p.length+" fields in "+file); }
				if(!("native".equals(p[0]) || "foreign".equals(p[0]))) {
					throw new RuntimeException("Invalid source-manifest role in "+file+": "+p[0]);
				}
				for(String value : p) {
					if(value.length()==0 || value.indexOf('\r')>=0 || value.indexOf('\n')>=0) {
						throw new RuntimeException("Empty/newline source-manifest field in "+file);
					}
				}
				final SourceRow row=new SourceRow(p[0],p[1],p[2],p[3],p[4]);
				if(!keys.add(row.key())) { throw new RuntimeException("Duplicate source role/assembly in "+file+": "+row.key()); }
				out.add(row);
			}
			if(!header || out.isEmpty()) { throw new RuntimeException("Source manifest has no rows: "+file); }
		} catch(IOException ex) { throw new RuntimeException("Could not read source manifest "+file,ex); }
		return out;
	}

	private static void validateSourceInputHashes(String sourceManifestPath, String sidecarPath,
			ArrayList<SourceRow> sources) {
		final HashMap<String,String[]> declared=new HashMap<String,String[]>();
		boolean sawHeader=false;
		try(BufferedReader br=new BufferedReader(new FileReader(sidecarPath))) {
			for(String line; (line=br.readLine())!=null;) {
				if(!line.startsWith("#source_input_sha256\t")) { continue; }
				final String[] p=line.split("\\t",-1);
				if(p.length==6 && "role".equals(p[1]) && "assembly_id".equals(p[2])) {
					if(sawHeader) { throw new RuntimeException("Duplicate source-input header in "+sidecarPath); }
					if(!line.equals("#source_input_sha256\t"+SOURCE_MANIFEST_COLUMNS)) {
						throw new RuntimeException("Invalid source-input header in "+sidecarPath);
					}
					sawHeader=true; continue;
				}
				if(p.length!=6 || !("native".equals(p[1]) || "foreign".equals(p[1]))
						|| !p[3].matches(HEX) || !p[4].matches(HEX) || !p[5].matches(HEX)) {
					throw new RuntimeException("Invalid source-input hash row in "+sidecarPath);
				}
				final String key=p[1]+"\t"+p[2];
				if(declared.put(key,new String[]{p[3],p[4],p[5]})!=null) {
					throw new RuntimeException("Duplicate source-input hash row in "+sidecarPath+": "+key);
				}
			}
		} catch(IOException ex) { throw new RuntimeException("Could not read sidecar source-input hashes "+sidecarPath,ex); }
		if(!sawHeader) { throw new RuntimeException("Sidecar lacks source-input hash header: "+sidecarPath); }
		final HashSet<String> expectedKeys=new HashSet<String>();
		for(SourceRow source : sources) {
			final String key=source.key();
			expectedKeys.add(key);
			final String[] hashes=declared.get(key);
			if(hashes==null) { throw new RuntimeException("Sidecar lacks source-input hash row: "+key); }
			final File gff=resolveExistingPath(sourceManifestPath,source.referenceGff,"reference GFF");
			final File assembly=resolveExistingPath(sourceManifestPath,source.assemblyFasta,"assembly FASTA");
			final File selected=resolveExistingPath(sourceManifestPath,source.selectedShreds,"selected shreds");
			if(!hashes[0].equalsIgnoreCase(MagQCBinManifest.sha256File(gff.getPath()))
					|| !hashes[1].equalsIgnoreCase(MagQCBinManifest.sha256File(assembly.getPath()))
					|| !hashes[2].equalsIgnoreCase(MagQCBinManifest.sha256File(selected.getPath()))) {
				throw new RuntimeException("Sidecar source-input hash mismatch for "+key);
			}
		}
		if(!expectedKeys.equals(declared.keySet())) {
			throw new RuntimeException("Sidecar/source-manifest source-input role set mismatch: "+sidecarPath);
		}
	}

	private static ArrayList<String> readSelectedHeaders(File file) {
		final ArrayList<String> headers=new ArrayList<String>();
		final HashSet<String> seen=new HashSet<String>();
		try(BufferedReader br=new BufferedReader(new FileReader(file))) {
			String current=null;
			for(String line; (line=br.readLine())!=null;) {
				if(line.startsWith(">")) {
					if(current!=null) { addHeader(headers,seen,current,file); }
					current=normalizeHeader(line.substring(1));
					if(current.length()==0) { throw new RuntimeException("Empty selected FASTA header: "+file); }
				} else if(current==null && line.trim().length()!=0) {
					throw new RuntimeException("Selected FASTA sequence precedes header: "+file);
				}
			}
			if(current!=null) { addHeader(headers,seen,current,file); }
		} catch(IOException ex) { throw new RuntimeException("Could not read selected FASTA "+file,ex); }
		return headers;
	}

	private static void addHeader(ArrayList<String> headers, HashSet<String> seen, String header, File file) {
		if(!seen.add(header)) { throw new RuntimeException("Duplicate selected FASTA header: "+header+" in "+file); }
		headers.add(header);
	}

	private static File resolveExistingPath(String manifestPath, String raw, String label) {
		final File direct=new File(raw);
		if(direct.isFile()) { return direct; }
		final File parent=new File(manifestPath).getAbsoluteFile().getParentFile();
		final File relative=new File(parent,raw);
		if(relative.isFile()) { return relative; }
		throw new RuntimeException("Missing "+label+" "+raw+" from source manifest "+manifestPath);
	}

	private static String normalizeHeader(String raw) {
		final int tab=raw.indexOf('\t');
		final String first=tab>=0 ? raw.substring(0,tab) : raw;
		return first.replace(' ','_').trim();
	}

	private static void requireSafeSpanHeader(String header, boolean encoded) {
		if(header.length()==0 || (!encoded && header.indexOf(';')>=0) || header.indexOf('\t')>=0
				|| header.indexOf('\n')>=0 || header.indexOf('\r')>=0) {
			throw new RuntimeException("Unsafe selected span/header: "+header);
		}
	}

	private static String join(Iterable<String> values) {
		final StringBuilder out=new StringBuilder();
		for(String value : values) {
			if(out.length()>0) { out.append(';'); }
			out.append(value);
		}
		return out.length()==0 ? "-" : out.toString();
	}

	private static String sidecarSourceManifestHash(String sidecar) {
		String found=null;
		try(BufferedReader br=new BufferedReader(new FileReader(sidecar))) {
			for(String line; (line=br.readLine())!=null;) {
				if(line.startsWith("#source_manifest_sha256\t")) {
					if(found!=null) { throw new RuntimeException("Duplicate #source_manifest_sha256 in "+sidecar); }
					final String[] p=line.split("\\t",-1);
					if(p.length!=2 || p[1].length()==0) { throw new RuntimeException("Malformed source-manifest hash in "+sidecar); }
					found=p[1];
				}
			}
		} catch(IOException ex) { throw new RuntimeException("Could not read sidecar "+sidecar,ex); }
		if(found==null) { throw new RuntimeException("Sidecar lacks #source_manifest_sha256: "+sidecar); }
		return found;
	}

	/** Canonical v1 binding hash shared with the opt-in consumer. */
	public static String selectionBindingSha256(String nativeIds, String foreignIds) {
		validateBindingField(nativeIds,"native_source_ids");
		validateBindingField(foreignIds,"foreign_source_ids");
		return sha256Text("selection_binding_v1\nnative\t"+nativeIds+"\nforeign\t"+foreignIds+"\n");
	}

	/**
	 * Re-derives and checks the role-separated selected FASTA header set for a
	 * consumer replay.  This is deliberately public so an index loaded earlier
	 * cannot become stale merely because a selected FASTA was edited afterward.
	 */
	public static void validateSelectedSourceSpanIds(String sourceManifestPath, String sidecarPath,
			MagQCBinManifest.Bin bin, String expected) {
		if(expected==null || expected.length()==0) {
			throw new IllegalArgumentException("Missing expected selected source span IDs for "+bin.id);
		}
		final String actual=deriveSelectedSourceSpanIds(sourceManifestPath, sidecarPath, bin);
		if(!actual.equals(expected)) {
			throw new IllegalArgumentException("Selected source span/header set mismatch for "+bin.id);
		}
	}

	private static void validateBindingField(String value, String field) {
		if(value==null || value.length()==0 || value.indexOf('\t')>=0 || value.indexOf('\n')>=0 || value.indexOf('\r')>=0) {
			throw new RuntimeException("Invalid "+field+" binding");
		}
	}

	private static String sha256Text(String text) {
		try {
			final MessageDigest md=MessageDigest.getInstance("SHA-256");
			final byte[] digest=md.digest(text.getBytes(StandardCharsets.UTF_8));
			final StringBuilder sb=new StringBuilder(64);
			for(byte b : digest) { sb.append(String.format("%02x",b)); }
			return sb.toString();
		} catch(NoSuchAlgorithmException ex) { throw new RuntimeException("SHA-256 unavailable",ex); }
	}

	private static String hashExisting(String file, String label) {
		if(!new File(file).isFile()) { throw new RuntimeException("Missing "+label+": "+file); }
		return MagQCBinManifest.sha256File(file);
	}

	private static Map<String,String> parseArgs(String[] args, int start) {
		final HashMap<String,String> out=new HashMap<String,String>();
		for(int i=start;i<args.length;i++) {
			final int eq=args[i].indexOf('=');
			if(eq<=0 || eq==args[i].length()-1) { throw new RuntimeException("Expected key=value, got "+args[i]); }
			out.put(args[i].substring(0,eq),args[i].substring(eq+1));
		}
		return out;
	}

	private static String require(Map<String,String> p, String key) {
		final String value=p.get(key);
		if(value==null || value.length()==0) { throw new RuntimeException("Missing "+key+" argument"); }
		return value;
	}

	private static void fail(String file, String message) { throw new RuntimeException("Invalid sidecar index "+file+": "+message); }
}
