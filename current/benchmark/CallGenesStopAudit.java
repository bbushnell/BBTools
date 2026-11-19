package benchmark;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import gff.GffLine;
import parse.LineParser1;
import parse.Parse;
import parse.PreParser;
import prok.ProkObject;
import prot.DigestSuffix;
import shared.Shared;
import structures.ByteBuilder;
import structures.StringNum;

/** Retained-GFF audit; never classifies or calls genes. Reference matching follows
 * CompareGff.processInner/processLine, including last-wins collisions and flags.
 * @author Sayu
 */
public final class CallGenesStopAudit{

	/** Uses the caller's CURRENT type flags; does not change them. Coordinates are
	 * GffLine one-based trueStop, not the scanner's zero-based stop-codon left edge. */
	public static ReferenceSet loadReference(String gff){
		require(gff!=null, "Reference GFF path is required");
		return fromLines(GffLine.loadGffFile(gff, "CDS,rRNA,tRNA", true));
	}

	/** Input is already filtered exactly as GffLine.loadGffFile(...,true) filters it. */
	public static ReferenceSet fromLines(List<GffLine> lines){
		require(lines!=null, "Reference features are required");
		return new ReferenceSet(lines);
	}

	/** Read-only view of retained references; feature objects must not be mutated. */
	public static final class ReferenceSet{
		private ReferenceSet(List<GffLine> lines){
			sourceRows=lines.size();
			int duplicates=0;
			source=Collections.unmodifiableList(new ArrayList<GffLine>(lines));
			for(GffLine g:lines){
				require(g!=null && g.seqid!=null && g.type!=null, "Reference key fields must exist");
				if(map.put(new StringNum(g.seqid, g.trueStop()), g)!=null){duplicates++;}
			}
			duplicateKeys=duplicates;
			all=Collections.unmodifiableList(new ArrayList<GffLine>(map.values()));
			final HashMap<String, ArrayList<GffLine>> groups=new HashMap<String, ArrayList<GffLine>>();
			for(GffLine g:all){
				ArrayList<GffLine> a=groups.get(g.seqid);
				if(a==null){
					a=new ArrayList<GffLine>();
					groups.put(g.seqid, a);
				}
				a.add(g);
			}
			for(Map.Entry<String, ArrayList<GffLine>> e:groups.entrySet()){
				byContig.put(e.getKey(), Collections.unmodifiableList(e.getValue()));
			}
		}
		/** Return the last reference row for a contig/stop key, or null if absent. */
		public GffLine lookup(String seqid, int trueStop){return map.get(new StringNum(seqid, trueStop));}
		/** Test whether a query row matches the retained reference type, strand and true stop. */
		public boolean matches(GffLine query){
			require(query!=null, "Query feature is required");
			final GffLine ref=lookup(query.seqid, query.trueStop());
			return ref!=null && ref.strand==query.strand && ref.type.equals(query.type);
		}
		/** Return the number of unique contig/stop keys after last-wins duplicate handling. */
		public int uniqueStops(){return map.size();}
		/** Return every input reference row, including duplicate keys. */
		public List<GffLine> sourceRecords(){return source;}
		/** Return last-wins reference rows in stable insertion order. */
		public List<GffLine> records(){return all;}
		/** Return last-wins reference rows for one contig. */
		public List<GffLine> records(String seqid){
			final List<GffLine> a=byContig.get(seqid);
			return a==null ? Collections.<GffLine>emptyList() : a;
		}
		public final int sourceRows, duplicateKeys;
		private final LinkedHashMap<StringNum, GffLine> map=new LinkedHashMap<StringNum, GffLine>();
		private final HashMap<String, List<GffLine>> byContig=new HashMap<String, List<GffLine>>();
		private final List<GffLine> all, source;
	}

	/** Same eight reported integer counts as CompareGff: TPstart/TPstop/FNstart/
	 * FNstop reference-relative, then TPstart/TPstop/FPstart/FPstop query-relative.
	 * TP values count matching query rows; FN values count unseen reference keys. */
	public static long[] grade(ReferenceSet refs, List<GffLine> queries){
		require(refs!=null && queries!=null, "Reference and query features are required");
		final long[] n=new long[8];
		final HashSet<StringNum> starts=new HashSet<StringNum>(), stops=new HashSet<StringNum>();
		for(GffLine q:queries){
			if(!ProkObject.processType(q.prokType())){continue;}
			if(!refs.matches(q)){
				n[6]++;
				n[7]++;
				continue;
			}
			final StringNum key=new StringNum(q.seqid, q.trueStop());
			stops.add(key);
			n[1]++;
			n[5]++;
			if(refs.lookup(q.seqid, q.trueStop()).trueStart()==q.trueStart()){
				starts.add(key);
				n[0]++;
				n[4]++;
			}else{n[6]++;}
		}
		n[2]=refs.uniqueStops()-starts.size();
		n[3]=refs.uniqueStops()-stops.size();
		return n;
	}

	/** Exact semicolon attribute; the initial native comma-separated description is ignored. */
	public static String attribute(GffLine g, String key){
		require(g!=null && key!=null && !key.isEmpty(), "Feature and attribute key are required");
		if(g.attributes==null){return null;}
		final String prefix=key+"=";
		String value=null;
		for(int start=0;start<g.attributes.length();){
			int end=g.attributes.indexOf(';', start);
			if(end<0){end=g.attributes.length();}
			if(g.attributes.startsWith(prefix, start)){
				require(value==null, "Duplicate GFF attribute "+key);
				value=g.attributes.substring(start+prefix.length(), end);
			}
			start=end+1;
		}
		return value;
	}

	/** Audit a completed retained benchmark root and emit grade, support and length tables. */
	public static void main(String[] args) throws Exception{
		final PreParser pp=new PreParser(args, CallGenesStopAudit.class, false);
		try{
			Path root=null, out=null;
			for(String arg:pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				require(pair.length==2, "Require key=value");
				if(pair[0].equals("root")){
					require(root==null, "Duplicate root");
					root=Paths.get(pair[1]);
				}else if(pair[0].equals("out")){
					require(out==null, "Duplicate out");
					out=Paths.get(pair[1]);
				}else{throw new IllegalArgumentException("Unknown option "+pair[0]);}
			}
			require(root!=null && out!=null && !Files.exists(out), "Require completed root and fresh audit output");
			// Reproduce the benchmark's cdsonly bootstrap in this fresh, audit-only JVM.
			ProkObject.callCDS=true;
			ProkObject.calltRNA=false;
			ProkObject.call16S=false;
			ProkObject.call23S=false;
			ProkObject.call5S=false;
			ProkObject.call18S=false;
			run(root, out);
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(Path root, Path out) throws Exception{
		final HashMap<String, long[]> expected=readGrades(root.resolve("collection_all/grades.tsv"));
		require(expected.size()==900, "Retained benchmark must contain 900 unique grade rows");
		Files.createDirectory(out);
		final ByteStreamWriter summary=writer(out.resolve("selected_audit.tsv")),
			missed=writer(out.resolve("native_missed.tsv")), lengths=writer(out.resolve("annotation_lengths.tsv"));
		summary.println("genome\tarm\tref_rows\tref_keys\tduplicate_ref_keys\ttp_start\ttp_stop\tfn_start\tfn_stop\tfp_start\tfp_stop\tselected_supported\tsupported_annotated\tsupported_unannotated\tnative_score_n\tnative_score_min\tnative_score_median\tnative_score_max\tverified_bonus");
		missed.println("genome\tseqid\ttrue_stop_1based\tstrand\ttype\tstart_1based\tend_1based\tpartial\texception\tcalled400\tcalled400_supported\tunselected_scan_accept_candidate_evidence");
		lengths.println("genome\tscope\tCDS_rows\tCDS_lt600nt\tCDS_lt603nt\tCDS_partial\tCDS_partial_lt603nt\tCDS_length_not_multiple_of3");
		final HashSet<String> seen=new HashSet<String>();
		long nativeFN=0, nativeFP=0;
		try{
			for(int batch=0;batch<4;batch++){
				final Path job=root.resolve(Integer.toString(batch));
				pin(job.resolve("genomes.tsv"),
					new String(Files.readAllBytes(job.resolve("genomes.sha80")), StandardCharsets.US_ASCII).trim());
				final HashMap<String, String> packed=readPins(job.resolve("packaged.tsv"));
				final ByteFile in=ByteFile.makeByteFile1(job.resolve("genomes.tsv").toString(), false);
				final LineParser1 lp=new LineParser1('\t');
				try{
					byte[] line=in.nextLine();
					require(line!=null, "Genome manifest header missing");
					for(line=in.nextLine();line!=null;line=in.nextLine()){
						lp.set(line);
						require(lp.terms()==5, "Genome manifest requires five columns");
						final String genome=lp.parseString(0), gff=lp.parseString(3);
						require(genome.matches("g[0-9]{4}") && seen.add(genome), "Unsafe or repeated genome key");
						pin(Paths.get(gff), lp.parseString(4));
						final ReferenceSet refs=loadReference(gff);
						final ArrayList<GffLine> nativeCalls=readCalls(job, genome, "native", packed),
							zero=readCalls(job, genome, "zero_effect", packed),
							bonus=readCalls(job, genome, "bonus400", packed);
						final long[] n=writeArm(summary, genome, "native", refs, nativeCalls, expected),
							z=writeArm(summary, genome, "zero_effect", refs, zero, expected);
						writeArm(summary, genome, "bonus400", refs, bonus, expected);
						require(Arrays.equals(n, z), "Retained zero-effect grade changed");
						nativeFN+=n[3];
						nativeFP+=n[7];
						final HashSet<StringNum> nativeStops=matchedKeys(refs, nativeCalls),
							bonusStops=matchedKeys(refs, bonus), supportedStops=new HashSet<StringNum>();
						writeLengths(lengths, genome, "source_rows", refs.sourceRecords(), null);
						writeLengths(lengths, genome, "grader_keys", refs.records(), null);
						writeLengths(lengths, genome, "native_missed_grader_keys", refs.records(), nativeStops);
						for(GffLine q:bonus){
							if(refs.matches(q) && supported(q)){supportedStops.add(new StringNum(q.seqid, q.trueStop()));}
						}
						for(GffLine ref:refs.records()){
							final StringNum key=new StringNum(ref.seqid, ref.trueStop());
							if(nativeStops.contains(key)){continue;}
							final ByteBuilder b=new ByteBuilder().append(genome).tab().append(ref.seqid)
								.tab().append(ref.trueStop()).tab().append(ref.strand).tab().append(ref.type)
								.tab().append(ref.start).tab().append(ref.stop).tab().append(value(attribute(ref, "partial")))
								.tab().append(value(attribute(ref, "exception"))).tab().append(bonusStops.contains(key) ? 1 : 0)
								.tab().append(supportedStops.contains(key) ? 1 : 0).tab().append("UNRECORDED").nl();
							missed.print(b);
						}
					}
				}finally{require(!in.close(), "Genome manifest read failed");}
			}
			require(seen.size()==100 && nativeFN==23008 && nativeFP==12172,
				"Retained native reconciliation failed: genomes="+seen.size()+" FN="+nativeFN+" FP="+nativeFP);
		}finally{
			final boolean bad=summary.poisonAndWait() | missed.poisonAndWait() | lengths.poisonAndWait();
			require(!bad, "Audit output write failed");
		}
		Files.write(out.resolve("PASS"),
			"CALLGENES_STOP_AUDIT_PASS genomes=100 native_fn=23008 native_fp=12172\n"
			.getBytes(StandardCharsets.US_ASCII));
		System.out.println("CALLGENES_STOP_AUDIT_PASS genomes=100 native_fn=23008 native_fp=12172");
	}
	/** Annotation lengths only: scanner eligibility depends on uninterrupted upstream
	 * codons, not annotated start or partial flags. Emit both 600nt and 603nt cutoffs;
	 * a complete CDS normally includes the 3nt stop, whereas a partial may not. */
	private static void writeLengths(ByteStreamWriter out, String genome, String scope,
		List<GffLine> records, HashSet<StringNum> exclude){
		final long[] n=new long[6];
		for(GffLine g:records){
			if(!g.type.equals("CDS") || (exclude!=null && exclude.contains(new StringNum(g.seqid, g.trueStop())))){continue;}
			require(g.length()>0, "Reference CDS length must be positive");
			n[0]++;
			if(g.length()<600){n[1]++;}
			if(g.length()<603){n[2]++;}
			if(g.partial()){
				n[3]++;
				if(g.length()<603){n[4]++;}
			}
			if(g.length()%3!=0){n[5]++;}
		}
		final ByteBuilder b=new ByteBuilder().append(genome).tab().append(scope);
		for(long v:n){b.tab().append(v);}
		out.print(b.nl());
	}

	private static long[] writeArm(ByteStreamWriter out, String genome, String arm, ReferenceSet refs,
		List<GffLine> calls, HashMap<String, long[]> expected){
		final long[] n=grade(refs, calls);
		require(Arrays.equals(n, expected.get(genome+"/"+arm)),
			"Eight-count grader mismatch "+genome+"/"+arm+": "+Arrays.toString(n));
		long support=0, annotated=0;
		final double[] scores=new double[calls.size()];
		int count=0;
		final int bonus=arm.equals("bonus400") ? 400 : 0;
		for(GffLine q:calls){
			if(supported(q)){
				final String s=attribute(q, "bbgap_native_score"), b=attribute(q, "bbgap_call_bonus");
				require(s!=null && b!=null, "Supported call lacks score/bonus provenance");
				final double score=Double.parseDouble(s), actual=Double.parseDouble(b);
				require(Double.isFinite(score) && actual==bonus,
					"Nonfinite native score or wrong applied bonus "+genome+"/"+arm);
				scores[count++]=score;
				support++;
				if(refs.matches(q)){annotated++;}
			}
		}
		Arrays.sort(scores, 0, count);
		final ByteBuilder b=new ByteBuilder().append(genome).tab().append(arm).tab().append(refs.sourceRows)
			.tab().append(refs.uniqueStops()).tab().append(refs.duplicateKeys);
		for(int slot:new int[]{0, 1, 2, 3, 6, 7}){b.tab().append(n[slot]);}
		b.tab().append(support).tab().append(annotated).tab().append(support-annotated).tab().append(count);
		if(count==0){
			b.append("\tNA\tNA\tNA\tNA");
		}else{
			b.tab().appendSlow(scores[0]).tab().appendSlow((scores[(count-1)/2]+scores[count/2])*0.5)
				.tab().appendSlow(scores[count-1]).tab().append(bonus);
		}
		out.print(b.nl());
		return n;
	}
	private static boolean supported(GffLine q){return "stop_supported".equals(attribute(q, "bbgap_status"));}
	private static String value(String s){return s==null?"absent":s;}
	private static HashSet<StringNum> matchedKeys(ReferenceSet refs, List<GffLine> calls){
		final HashSet<StringNum> keys=new HashSet<StringNum>();
		for(GffLine q:calls){
			if(ProkObject.processType(q.prokType()) && refs.matches(q)){keys.add(new StringNum(q.seqid, q.trueStop()));}
		}
		return keys;
	}
	private static ArrayList<GffLine> readCalls(Path job, String genome, String arm,
		HashMap<String, String> pins) throws Exception{
		final String relative=genome+"/"+arm+".gff.gz";
		final Path file=job.resolve("output").resolve(relative);
		pin(file, pins.get(relative));
		return GffLine.loadGffFile(file.toString(), null, false);
	}
	private static HashMap<String, long[]> readGrades(Path path){
		final HashMap<String, long[]> rows=new HashMap<String, long[]>();
		final ByteFile in=ByteFile.makeByteFile1(path.toString(), false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			require(in.nextLine()!=null, "Grade header missing");
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				lp.set(line);
				require(lp.terms()==10, "Expected ten grade columns");
				final long[] n=new long[8];
				for(int i=0;i<8;i++){n[i]=lp.parseLong(i+2);}
				require(rows.put(lp.parseString(0)+"/"+lp.parseString(1), n)==null, "Repeated grade key");
			}
		}finally{require(!in.close(), "Grade table read failed");}
		return rows;
	}
	private static HashMap<String, String> readPins(Path path){
		final HashMap<String, String> rows=new HashMap<String, String>();
		final ByteFile in=ByteFile.makeByteFile1(path.toString(), false);
		final LineParser1 lp=new LineParser1('\t');
		try{
			for(byte[] line=in.nextLine();line!=null;line=in.nextLine()){
				lp.set(line);
				require(lp.terms()==2 && rows.put(lp.parseString(0), lp.parseString(1))==null,
					"Invalid packaged receipt");
			}
		}finally{require(!in.close(), "Packaged receipt read failed");}
		return rows;
	}
	private static void pin(Path file, String expected) throws Exception{
		require(expected!=null && expected.matches("[0-9a-f]{20}") && Files.isRegularFile(file)
			&& !Files.isSymbolicLink(file) && expected.equals(DigestSuffix.file(file.toString())),
			"Input sha80 mismatch: "+file);
	}
	private static ByteStreamWriter writer(Path path){
		final ByteStreamWriter w=new ByteStreamWriter(path.toString(), false, false, false);
		w.start();
		return w;
	}
	private static void require(boolean ok, String message){if(!ok){throw new IllegalArgumentException(message);}}
}
