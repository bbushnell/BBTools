package prot;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/**
 * Creates the reference-CDS interval/survival sidecar used by the redesigned gene
 * labels.  It deliberately does not read fresh CallGenes CDS totals or family hits:
 * those remain observable features, while this file is reference-identity truth.
 *
 * <p>Input manifest columns are:
 * {@code role<TAB>assembly_id<TAB>reference_gff<TAB>assembly_fasta<TAB>selected_shreds}.
 * Role is {@code native} or {@code foreign}.  GFF3 CDS coordinates are 1-based
 * inclusive.  BBTools shred suffixes are parsed as 0-based inclusive start-stop and
 * converted internally to half-open intervals.  A gene survives only when one
 * selected source span contains its complete interval.</p>
 *
 * <p>The tool is intentionally project-only for this increment.  The class is NEW
 * under {@code src/prot/} and is suitable for Brian's later BBTools copy after review.</p>
 */
public final class ReferenceCdsSurvivalSidecar {

	/** v1 = raw identifiers (a ',' inside a shred id was unrepresentable in containing_shred_ids); v2 (2026-09-09) = every
	 *  containing_shred_ids item is IdentifierCodec-encoded (reserved bytes tab , ; | CR LF % as %XX), gene_key stays opaque. */
	static final String SCHEMA_V1="reference_cds_survival_v1", SCHEMA="reference_cds_survival_v2";
	// Package-private (UMP45 2026-09-09): shared verbatim with ReferenceCdsShredSurvivalTable so the
	// per-assembly precomputation parses shred ids and gene identities by exactly this code path.
	static final Pattern SHRED=Pattern.compile("^(.+)_([0-9]+)-([0-9]+)(?:_tid_[0-9]+)?$");
	/** v2 identifier components are IdentifierCodec-escaped at output, so the instance identity
	 * may contain delimiter and percent bytes before that encoding step. */
	static final Pattern INSTANCE_SUFFIX=Pattern.compile("^(.*)__instance_(.+)$");
	private static final String MANIFEST_HEADER="role\tassembly_id\treference_gff\tassembly_fasta\tselected_shreds";

	private ReferenceCdsSurvivalSidecar(){}

	private static final class Source {
		final String role, assemblyId, gff, assemblyFasta, selectedShreds;
		Source(String role_, String assemblyId_, String gff_, String assemblyFasta_, String selectedShreds_){
			role=role_; assemblyId=assemblyId_; gff=gff_; assemblyFasta=assemblyFasta_; selectedShreds=selectedShreds_;
		}
	}

	private static final class Span {
		final String source, header;
		final int start, end;
		Span(String source_, String header_, int start_, int end_){source=source_; header=header_; start=start_; end=end_;}
		String key(){return source+"\t"+start+"\t"+end;}
	}

	private static final class Gene {
		final Source source;
		final String contig, id, strand;
		int start=Integer.MAX_VALUE, end=-1;
		Gene(Source source_, String contig_, String id_, String strand_){source=source_; contig=contig_; id=id_; strand=strand_;}
		String key(){return source.role+"|"+source.assemblyId+"|"+contig+"|"+id+"|"+(start+1)+"|"+end+"|"+strand;}
	}

	private static final class SourceData {
		final Source source;
		final TreeMap<String,Integer> lengths=new TreeMap<String,Integer>();
		final ArrayList<Span> spans=new ArrayList<Span>();
		final TreeMap<String,Gene> genes=new TreeMap<String,Gene>();
		final ArrayList<String> hashes=new ArrayList<String>();
		SourceData(Source source_){source=source_;}
	}

	static final class Counts { int genes, retained; }

	public static void main(final String[] args) throws Exception{
		if(args.length==1 && args[0].equalsIgnoreCase("selftest")){
			System.err.println("ReferenceCdsSurvivalSidecar: use ReferenceCdsSurvivalSidecarTest for selftest");
			return;
		}
		String manifest=null, out=null;
		for(String arg : args){
			int eq=arg.indexOf('=');
			if(eq<1){throw new RuntimeException("Arguments must be key=value: "+arg);}
			String key=arg.substring(0, eq).toLowerCase(), value=arg.substring(eq+1);
			if(key.equals("manifest")){manifest=value;}
			else if(key.equals("out")){out=value;}
			else{throw new RuntimeException("Unknown argument: "+arg);}
		}
		if(manifest==null || out==null){throw new RuntimeException("Usage: java -ea prot.ReferenceCdsSurvivalSidecar manifest=<source_manifest.tsv> out=<sidecar.tsv>");}
		Counts counts=write(manifest, out);
		System.err.println("ReferenceCdsSurvivalSidecar PASS: genes="+counts.genes+" retained="+counts.retained+" out="+out);
	}

	static Counts write(final String manifestPath, final String outPath) throws Exception{
		final List<Source> sources=readManifest(manifestPath);
		final ArrayList<SourceData> data=new ArrayList<SourceData>();
		final HashSetLike sourceKeys=new HashSetLike();
		for(Source source : sources){
			if(!sourceKeys.add(source.role+"\t"+source.assemblyId)){
				throw new IllegalArgumentException("Duplicate role/assembly_id in manifest: "+source.role+"/"+source.assemblyId);
			}
			SourceData d=new SourceData(source);
			d.hashes.add(sha256(source.gff));
			d.hashes.add(sha256(source.assemblyFasta));
			d.hashes.add(sha256(source.selectedShreds));
			readAssembly(d);
			readSelectedShreds(d);
			readGff(d);
			data.add(d);
		}

		final Path outFile=Paths.get(outPath);
		Path parent=outFile.toAbsolutePath().getParent();
		if(parent!=null){Files.createDirectories(parent);}
		final Counts total=new Counts();
		try(BufferedWriter w=Files.newBufferedWriter(outFile, StandardCharsets.UTF_8)){
			w.write("#schema_version\t"+SCHEMA); w.newLine();
			w.write("#tool\tprot.ReferenceCdsSurvivalSidecar"); w.newLine();
			w.write("#coordinate_convention\tGFF3=1-based-inclusive; shred_suffix=0-based-inclusive; internal=0-based-half-open"); w.newLine();
			w.write("#header_normalization\tfirst-token-before-tab; spaces-to-underscore; no other rewriting"); w.newLine();
			w.write("#gene_universe\tprotein-coding CDS rows from original reference GFF; RNA and fresh calls excluded"); w.newLine();
			w.write("#gene_id_attribute_priority\tgene_id,locus_tag,Parent,ID; multiple CDS rows with one identity form one interval envelope"); w.newLine();
			w.write("#gene_id_fallback\tattribute-less CDS uses its 1-based physical GFF line number as identity; multipart CDS cannot be inferred and remains one row per CDS feature"); w.newLine();
			w.write("#survival_rule\tfull reference interval contained by one selected source span; split pieces never survive"); w.newLine();
			w.write("#status_values\tRETAINED=one full containing span; SPLIT_OR_PARTIAL=selected overlap but no full span; OMITTED=no selected overlap"); w.newLine();
			w.write("#instance_identity\tselected shred header is the instance identity; optional suffix=__instance_<nonempty-instance-id>; v2 output percent-escapes reserved identifier bytes; repeated source spans are allowed only with distinct headers; exact duplicate headers fail"); w.newLine();
			w.write("#zero_reference_cds_allowed\ttrue; a source with zero CDS rows emits no gene rows and leaves the downstream 0/0 convention to its consumer"); w.newLine();
			w.write("#fresh_called_cds_counts\tfeature_only_not_a_label_operand"); w.newLine();
			w.write("#mapped_family_counts\tfeature_only_not_a_label_operand"); w.newLine();
			w.write("#source_manifest_sha256\t"+sha256(manifestPath)); w.newLine();
			w.write("#source_count\t"+data.size()); w.newLine();
			w.write("#source_input_sha256\trole\tassembly_id\treference_gff\tassembly_fasta\tselected_shreds"); w.newLine();
			for(SourceData d : data){
				w.write("#source_input_sha256\t"+d.source.role+"\t"+d.source.assemblyId+"\t"+d.hashes.get(0)+"\t"+d.hashes.get(1)+"\t"+d.hashes.get(2)); w.newLine();
			}
			w.write("role\tassembly_id\tgene_key\tsource_contig\tgene_id\tstart1\tend1\tstrand\tstatus\tcontaining_shred_count\tcontaining_shred_ids"); w.newLine();
			for(SourceData d : data){
				for(Gene gene : d.genes.values()){
					total.genes++;
					ArrayList<String> containing=new ArrayList<String>();
					for(Span span : d.spans){
						if(span.source.equals(gene.contig) && span.start<=gene.start && span.end>=gene.end){containing.add(IdentifierCodec.encode(span.header));}
					}
					Collections.sort(containing);
					if(!containing.isEmpty()){total.retained++;}
					boolean overlap=false;
					for(Span span : d.spans){
						if(span.source.equals(gene.contig) && span.start<gene.end && span.end>gene.start){overlap=true; break;}
					}
					String status=(containing.isEmpty() ? (overlap ? "SPLIT_OR_PARTIAL" : "OMITTED") : "RETAINED");
					w.write(gene.source.role+"\t"+gene.source.assemblyId+"\t"+gene.key()+"\t"+gene.contig+"\t"+gene.id+"\t"+(gene.start+1)+"\t"+gene.end+"\t"+gene.strand+"\t"+status+"\t"+containing.size()+"\t"+join(containing));
					w.newLine();
				}
			}
		}
		return total;
	}

	private static List<Source> readManifest(final String path) throws IOException{
	ArrayList<Source> list=new ArrayList<Source>();
		try(BufferedReader r=reader(path)){
			String line; boolean header=false; int lineNo=0;
			while((line=r.readLine())!=null){
				lineNo++;
				if(line.isEmpty() || line.startsWith("#")){continue;}
				if(!header){
					if(!line.equals(MANIFEST_HEADER)){throw new IllegalArgumentException("Manifest header must be exactly: "+MANIFEST_HEADER);}
					header=true; continue;
				}
				String[] f=line.split("\t", -1);
				if(f.length!=5){throw new IllegalArgumentException("Manifest line "+lineNo+" has "+f.length+" fields, expected 5");}
				if(!(f[0].equals("native") || f[0].equals("foreign"))){throw new IllegalArgumentException("Manifest line "+lineNo+" role must be native or foreign");}
				for(int i=1;i<5;i++){if(f[i].isEmpty()){throw new IllegalArgumentException("Manifest line "+lineNo+" has empty field "+i);}}
				list.add(new Source(f[0], f[1], f[2], f[3], f[4]));
			}
			if(!header || list.isEmpty()){throw new IllegalArgumentException("Manifest has no source rows: "+path);}
		}
		return list;
	}

	private static void readAssembly(final SourceData d) throws IOException{
		try(BufferedReader r=reader(d.source.assemblyFasta)){
			String line, key=null; int len=0; int lineNo=0;
			while((line=r.readLine())!=null){
				lineNo++;
				if(line.startsWith(">")){
					if(key!=null){putLength(d, key, len);}
					key=normalizeHeader(line.substring(1)); len=0;
					if(key.isEmpty()){throw new IllegalArgumentException("Empty assembly FASTA header at "+d.source.assemblyFasta+":"+lineNo);}
				}else if(key!=null){len+=sequenceLength(line);}
				else if(!line.trim().isEmpty()){throw new IllegalArgumentException("Sequence before header at "+d.source.assemblyFasta+":"+lineNo);}
			}
			if(key!=null){putLength(d, key, len);}
		}
		if(d.lengths.isEmpty()){throw new IllegalArgumentException("Assembly FASTA has no records: "+d.source.assemblyFasta);}
	}

	private static void readSelectedShreds(final SourceData d) throws IOException{
		final Map<String,Boolean> seenHeaders=new TreeMap<String,Boolean>();
		try(BufferedReader r=reader(d.source.selectedShreds)){
			String line, header=null, source=null; int start=0, end=0, len=0, lineNo=0;
			while((line=r.readLine())!=null){
				lineNo++;
				if(line.startsWith(">")){
					if(header!=null){addSpan(d, seenHeaders, header, source, start, end, len, lineNo);}
					header=normalizeHeader(line.substring(1));
					String coordinateHeader=header;
					Matcher instance=INSTANCE_SUFFIX.matcher(coordinateHeader);
					if(instance.matches()){coordinateHeader=instance.group(1);}
					Matcher m=SHRED.matcher(coordinateHeader);
					if(!m.matches()){throw new IllegalArgumentException("Cannot parse shred source span header at "+d.source.selectedShreds+":"+lineNo+": "+header);}
					source=normalizeHeader(m.group(1)); start=parseInt(m.group(2), "shred start");
					int stop=parseInt(m.group(3), "shred stop");
					if(stop<start){throw new IllegalArgumentException("Shred stop before start: "+header);}
					end=stop+1; len=0;
				}else if(header!=null){len+=sequenceLength(line);}
				else if(!line.trim().isEmpty()){throw new IllegalArgumentException("Sequence before shred header at "+d.source.selectedShreds+":"+lineNo);}
			}
			if(header!=null){addSpan(d, seenHeaders, header, source, start, end, len, lineNo);}
		}
	}

	private static void addSpan(SourceData d, Map<String,Boolean> seenHeaders,
			String header, String source, int start, int end, int sequenceLength, int lineNo){
		if(seenHeaders.put(header, Boolean.TRUE)!=null){throw new IllegalArgumentException("Duplicate selected shred header: "+header);}
		if(sequenceLength!=end-start){throw new IllegalArgumentException("Selected shred sequence length "+sequenceLength+" != encoded span length "+(end-start)+" at line "+lineNo+": "+header);}
		Integer length=d.lengths.get(source);
		if(length==null){throw new IllegalArgumentException("Selected shred source key not found in assembly FASTA: "+source+" (header "+header+")");}
		if(start<0 || end>length){throw new IllegalArgumentException("Selected shred span outside source contig "+source+": "+start+"-"+(end-1)+" length="+length);}
		d.spans.add(new Span(source, header, start, end));
	}

	private static void readGff(final SourceData d) throws IOException{
		try(BufferedReader r=reader(d.source.gff)){
			String line; int lineNo=0;
			while((line=r.readLine())!=null){
				lineNo++;
				if(line.isEmpty() || line.charAt(0)=='#'){continue;}
				String[] f=line.split("\t", -1);
				if(f.length<9){throw new IllegalArgumentException("GFF line "+lineNo+" has fewer than 9 fields: "+d.source.gff);}
				if(!f[2].equals("CDS")){continue;}
				String contig=normalizeHeader(f[0]);
				Integer length=d.lengths.get(contig);
				if(length==null){throw new IllegalArgumentException("GFF CDS seqid absent from assembly FASTA: "+contig+" line "+lineNo);}
				int start=parseInt(f[3], "GFF start"); int end=parseInt(f[4], "GFF end");
				if(start<1 || end<start || end>length){throw new IllegalArgumentException("GFF CDS interval outside assembly at line "+lineNo+": "+start+"-"+end+" length="+length);}
				String strand=f[6];
				if(!(strand.equals("+") || strand.equals("-") || strand.equals("."))){throw new IllegalArgumentException("Unsupported GFF strand at line "+lineNo+": "+strand);}
				String id=geneId(f[8]);
				if(id==null){id="CDS_row_"+lineNo;}
				String key=d.source.role+"|"+d.source.assemblyId+"|"+contig+"|"+id+"|"+strand;
				Gene gene=d.genes.get(key);
				if(gene==null){gene=new Gene(d.source, contig, id, strand); d.genes.put(key, gene);}
				gene.start=Math.min(gene.start, start-1); gene.end=Math.max(gene.end, end);
			}
		}
	}

	static String geneId(final String attrs){
		Map<String,String> map=new LinkedHashMap<String,String>();
		for(String item : attrs.split(";")){
			int eq=item.indexOf('='); if(eq<0){eq=item.indexOf(' ');}
			if(eq>0){String key=item.substring(0,eq).trim(); String value=item.substring(eq+1).trim(); if(!value.isEmpty()){map.put(key, value.replace(' ', '_'));}}
		}
		for(String key : new String[]{"gene_id","locus_tag","Parent","ID"}){String value=map.get(key); if(value!=null){return value.split(",",2)[0];}}
		return null;
	}

	private static void putLength(SourceData d, String key, int len){
		if(len<=0){throw new IllegalArgumentException("Empty assembly sequence for "+key);}
		if(d.lengths.put(key, len)!=null){throw new IllegalArgumentException("Duplicate assembly FASTA header: "+key);}
	}

	private static int sequenceLength(String line){
		int n=0; for(int i=0;i<line.length();i++){if(!Character.isWhitespace(line.charAt(i))){n++;}}
		return n;
	}

	static String normalizeHeader(String raw){
		String s=raw; int tab=s.indexOf('\t'); if(tab>=0){s=s.substring(0,tab);}
		return s.replace(' ', '_').trim();
	}

	private static int parseInt(String value, String what){
		try{return Integer.parseInt(value);}catch(NumberFormatException e){throw new IllegalArgumentException("Invalid "+what+": "+value, e);}
	}

	private static String join(List<String> values){
		if(values.isEmpty()){return "-";}
		StringBuilder b=new StringBuilder(); for(String s : values){if(b.length()>0){b.append(',');} b.append(s);}
		return b.toString();
	}

	private static BufferedReader reader(String path) throws IOException{
		InputStream in=new FileInputStream(path);
		if(path.toLowerCase().endsWith(".gz")){in=new GZIPInputStream(in, 1<<16);}
		return new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8), 1<<16);
	}

	static String sha256(String path) throws IOException{
		try{
			MessageDigest md=MessageDigest.getInstance("SHA-256");
			try(InputStream in=new FileInputStream(path)){
				byte[] buf=new byte[1<<16]; int n;
				while((n=in.read(buf))>=0){if(n>0){md.update(buf,0,n);}}
			}
			return hex(md.digest());
		}catch(java.security.NoSuchAlgorithmException e){throw new AssertionError(e);}
	}

	private static String hex(byte[] bytes){
		StringBuilder b=new StringBuilder(bytes.length*2); for(byte x : bytes){b.append(String.format("%02x", x & 0xff));} return b.toString();
	}

	/** Tiny local set with no dependency on project collection classes. */
	private static final class HashSetLike {
		private final java.util.HashSet<String> set=new java.util.HashSet<String>();
		boolean add(String s){return set.add(s);}
	}
}
