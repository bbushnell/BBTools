package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.Writer;
import stream.WriterFactory;
import structures.ByteBuilder;

/** Unions exact predicted loci for independent review without regrading any run.
 * Globally unique contig IDs bind each locus to one original genome. Overlapping
 * nonidentical calls and opposite strands remain separate, with all memberships.
 * @author Brian, Raiden
 */
public final class NcrnaReviewUnion {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaReviewUnion.class, false);
		String runs=null, roster=null, out=null, family=null;int genomes=-1, contigs=-1, count=-1, flank=-1;boolean panel=false;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("runs")){runs=value;}else if(key.equals("roster")){roster=value;}
			else if(key.equals("outprefix")){out=value;}else if(key.equals("family")){family=value;}
			else if(key.equals("genomes")){genomes=Integer.parseInt(value);}else if(key.equals("contigs")){contigs=Integer.parseInt(value);}
			else if(key.equals("runcount")){count=Integer.parseInt(value);}else if(key.equals("flank")){flank=Integer.parseInt(value);}
			else if(key.equals("panel")){panel=parse.Parse.parseBoolean(value);}
			else{throw new IllegalArgumentException("Unknown option: "+key);}
		}
		require(runs!=null && roster!=null && out!=null && family!=null && genomes>0 && contigs>0 && count>0 && flank>=0,
			"Bind reviewed run ledgers, original genome roster, family and extraction denominators");
		Shared.setThreads(1);Shared.TRIM_READ_DESCRIPTION=false;
		try{
			require(!panel || family.equals("cobalamin"), "Panel schema currently describes the three B12 subtypes only");
			final NcrnaReviewUnion x=new NcrnaReviewUnion(family, flank, panel);
			x.loadRuns(runs, count);x.extract(roster, out, genomes, contigs);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	private NcrnaReviewUnion(final String family_, final int flank_, final boolean panel_){family=family_;flank=flank_;panel=panel_;}

	/*--------------------------------------------------------------*/
	/*----------------            Inputs            ----------------*/
	/*--------------------------------------------------------------*/

	private void loadRuns(final String path, final int expected){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> names=new ObjectSet<String>(String.class), paths=new ObjectSet<String>(String.class);
		try{
			RrnaResourceIO.header(in, "run\tcalls\treview_candidates", path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==3, "Run manifest must retain its expected candidate count");
				final String run=p.parseString(0), calls=p.parseString(1);final int candidates=p.parseInt(2);
				require(!run.isEmpty() && names.add(run) && paths.add(calls) && candidates>=0, "Repeated run/path or invalid candidate count: "+run);
				if(panel){loadPanelCalls(run, calls, candidates);}else{loadCalls(run, calls, candidates);}
			}
		}finally{require(!in.close(), "Run manifest read failed");}
		require(names.size()==expected, "The union must consume every declared run exactly once");runCount=names.size();
		loci.sort((a, b)->{int n=a.seqid.compareTo(b.seqid);if(n==0){n=Character.compare(a.strand, b.strand);}
			if(n==0){n=Integer.compare(a.start, b.start);}return n==0?Integer.compare(a.stop, b.stop):n;});
		for(int i=0; i<loci.size(); i++){loci.get(i).id="review_"+i;}
	}
	private void loadCalls(final String run, final String path, final int expected){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);int retained=0;
		try{
			RrnaResourceIO.header(in, CALL_HEADER, path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==10, "Complete per-run disposition ledger required");
				final String call=p.parseString(0);require(ids.add(call), "Repeated candidate ID within run: "+run+" "+call);
				require(p.termEquals(family, 2) && (p.termEquals("true", 9) || p.termEquals("false", 9)), "Wrong family or malformed review flag");
				if(p.termEquals("false", 9)){continue;}
				require(p.termEquals("false", 8), "Exact duplicate query accounting stays in its original run, not the independent scan");
				final String seqid=p.parseString(1), strand=p.parseString(5);final int start=p.parseInt(3), stop=p.parseInt(4);
				require(start>0 && stop>=start && (strand.equals("+") || strand.equals("-")), "Review locus requires valid inclusive genomic coordinates");
				add(run, call, seqid, strand, start, stop, p.parseString(6), p.parseString(7), null);retained++;
			}
		}finally{require(!in.close(), "Candidate ledger read failed: "+path);}
		require(retained==expected, "Independent review candidates must match native grader accounting: "+run);
	}
	/** Selects only FULL_CM/union apparent FPs; preserves SOURCE outcomes without
	 * confusing subtype FPs (which include label swaps) with union FPs. */
	private void loadPanelCalls(final String run, final String path, final int expected){
		final ObjectMap<String,PanelCall> source=new ObjectMap<String,PanelCall>(String.class, PanelCall.class);
		final ObjectSet<String> full=new ObjectSet<String>(String.class);int retained=0;
		final LineParser1 p=new LineParser1('\t');
		for(int pass=0; pass<2; pass++){
			final ByteFile in=RrnaResourceIO.open(path);
			try{
				RrnaResourceIO.header(in, PANEL_CALL_HEADER, path);
				for(byte[] row; (row=in.nextLine())!=null; ){
					p.set(row);require(p.terms()==14, "Complete native panel audit schema required");
					if(!p.termEquals("union", 1) || !p.termEquals(pass==0?"SOURCE":"FULL_CM", 0)){continue;}
					final String id=p.parseString(2), outcome=p.parseString(10);final PanelCall c=new PanelCall(p);
					if(pass==0){
						require((outcome.equals("TP") || outcome.equals("APPARENT_FP")) && !source.contains(id), "SOURCE union call IDs/outcomes must be unique and scored");
						source.put(id, c);continue;
					}
					final PanelCall original=source.get(id);
					require(original!=null && original.same(c) && full.add(id), "SOURCE/FULL_CM must describe exactly the same call identities and coordinates: "+id);
					require(outcome.equals("TP") || outcome.equals("APPARENT_FP") || outcome.equals("UNSCORED"), "Unknown native panel outcome: "+outcome);
					if(!outcome.equals("APPARENT_FP")){continue;}
					add(run, id, c.seqid, c.strand, c.start, c.stop, original.outcome, outcome, c.genome);retained++;
				}
			}finally{require(!in.close(), "Panel audit read failed: "+path);}
		}
		require(source.size()==full.size() && retained==expected, "Every panel call needs both views and every union FP must enter independent review: "+run);
	}
	private void add(String run, String call, String seqid, String strand, int start, int stop, String source, String full, String genome){
		require(start>0 && stop>=start && (strand.equals("+") || strand.equals("-")), "Review locus requires valid inclusive genomic coordinates");
		final String key=seqid+'\t'+strand+'\t'+start+'\t'+stop;Locus l=physical.get(key);
		if(l==null){
			l=new Locus(seqid, strand.charAt(0), start, stop);l.expectedGenome=genome;physical.put(key, l);loci.add(l);
			Group g=byContig.get(seqid);if(g==null){g=new Group();byContig.put(seqid, g);}g.loci.add(l);
		}else{require(genome==null?l.expectedGenome==null:genome.equals(l.expectedGenome), "The same physical locus must retain its source genome across runs");}
		members.add(new Member(run, call, l, source, full));
	}

	/*--------------------------------------------------------------*/
	/*----------------          Extraction          ----------------*/
	/*--------------------------------------------------------------*/

	private void extract(final String roster, final String out, final int expectedGenomes, final int expectedContigs){
		RrnaResourceIO.fresh(out, ".contexts.fa.gz", ".loci.tsv", ".members.tsv", ".summary.tsv");
		final ObjectSet<String> genomes=new ObjectSet<String>(String.class), contigs=new ObjectSet<String>(String.class);
		final ByteFile in=RrnaResourceIO.open(roster);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, NcrnaCalibrationPrep.ROSTER_HEADER, roster);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==11, "Genome roster schema changed");
				final String genome=p.parseString(1), assembly=p.parseString(2), fasta=p.parseString(7);
				require(genomes.add(genome), "Duplicate source genome in extraction roster: "+genome);
				RrnaResourceIO.stream(fasta, r->{
					final String name=RrnaResourceIO.first(r.id);
					require(contigs.add(name), "Contig ID is ambiguous across original source genomes: "+name);
					final Group group=byContig.get(name);if(group==null){return;}
					require(r.bases!=null, "Review contig needs sequence: "+name);Tools.toUpperCase(r.bases);
					for(Locus l:group.loci){
						require(l.bases==null && l.stop<=r.length(), "Every review locus must bind exactly once inside its original contig");
						require(l.expectedGenome==null || genome.equals(l.expectedGenome), "Panel call must bind its original genome, not only a matching contig ID");
						final int first=(int)Math.max(0L, (long)l.start-1-flank), last=(int)Math.min(r.length()-1L, (long)l.stop-1+flank);
						final int left=l.start-1-first, right=last-l.stop+1;
						l.left=l.strand=='+'?left:right;l.right=l.strand=='+'?right:left;
						l.bases=Arrays.copyOfRange(r.bases, first, last+1);
						if(l.strand=='-'){AminoAcid.reverseComplementBasesInPlace(l.bases);}
						l.genome=genome;l.assembly=assembly;l.fasta=fasta;l.contigLength=r.length();
						assert(l.bases.length==l.stop-l.start+1+l.left+l.right):"CutGff-style clipped sense flanks must conserve context length";
					}
				});
			}
		}finally{require(!in.close(), "Genome roster read failed");}
		require(genomes.size()==expectedGenomes && contigs.size()==expectedContigs, "Frozen extraction genome/contig denominators changed");
		for(Locus l:loci){require(l.bases!=null, "Candidate contig absent from original genomes: "+l.seqid);}
		write(out, genomes.size(), contigs.size());
	}
	private void write(final String out, final int genomes, final int contigs){
		assert(genomes>0 && contigs>0):"Output is published only after original-source identity checks complete";
		final Writer fasta=WriterFactory.makeWriter(FileFormat.testOutput(out+".contexts.fa.gz", FileFormat.FASTA, null, true, false, false, true), null, 1);
		fasta.start();
		final ByteStreamWriter metadata=RrnaResourceIO.writer(out+".loci.tsv");
		final ByteStreamWriter membership=RrnaResourceIO.writer(out+".members.tsv"), summary=RrnaResourceIO.writer(out+".summary.tsv");
		final ByteBuilder b=new ByteBuilder();long bases=0;
		try{
			metadata.println(LOCUS_HEADER);membership.println("run\tcall_id\treview_id\tsource_outcome\tfull_outcome");
			final ArrayList<Read> records=new ArrayList<Read>(loci.size());
			for(Locus l:loci){
				records.add(new Read(l.bases, null, l.id+" lflank="+l.left+" rflank="+l.right, records.size(), false));
				metadata.print(b.clear().append(l.id).tab().append(l.genome).tab().append(l.assembly).tab().append(l.fasta).tab().append(l.seqid)
					.tab().append(family).tab().append(l.strand).tab().append(l.start).tab().append(l.stop).tab().append(l.contigLength)
					.tab().append(l.left).tab().append(l.right).tab().append(l.bases.length).append("\tPREDICTION_UNRESOLVED\n"));bases+=l.bases.length;
			}
			fasta.add(records, 0);
			for(Member m:members){membership.print(b.clear().append(m.run).tab().append(m.call).tab().append(m.locus.id).tab().append(m.source).tab().append(m.full).nl());}
			summary.println("runs\tgenomes\tcontigs\treview_events\tunique_loci\tcontext_bases\tflank");
			summary.print(b.clear().append(runCount).tab().append(genomes).tab().append(contigs).tab().append(members.size()).tab().append(loci.size()).tab().append(bases).tab().append(flank).nl());
		}finally{
			final boolean error=fasta.poisonAndWait();RrnaResourceIO.finish(metadata, membership, summary);
			require(!error, "Complete native FASTA output is required before reporting extraction success");
		}
		System.err.println("NCRNA_REVIEW_UNION_PASS events="+members.size()+" loci="+loci.size()+"; predictions only, no regrading or truth edits");
	}
	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Locus {
		Locus(String seqid_, char strand_, int start_, int stop_){seqid=seqid_;strand=strand_;start=start_;stop=stop_;}
		final String seqid;final char strand;final int start, stop;
		String id, genome, expectedGenome, assembly, fasta;int left, right, contigLength;byte[] bases;
	}
	private static final class PanelCall {
		PanelCall(LineParser1 p){
			genome=p.parseString(3);seqid=p.parseString(4);family=p.parseString(5);start=p.parseInt(6);stop=p.parseInt(7);strand=p.parseString(8);outcome=p.parseString(10);
			require(family.equals("RF00174") || family.equals("RF01482") || family.equals("RF01689"), "Unknown B12 subtype in panel audit");
		}
		boolean same(PanelCall c){return genome.equals(c.genome) && seqid.equals(c.seqid) && family.equals(c.family) && start==c.start && stop==c.stop && strand.equals(c.strand);}
		final String genome, seqid, family, strand, outcome;final int start, stop;
	}
	private static final class Group {final ArrayList<Locus> loci=new ArrayList<Locus>();}
	private static final class Member {
		Member(String run_, String call_, Locus locus_, String source_, String full_){run=run_;call=call_;locus=locus_;source=source_;full=full_;}
		final String run, call, source, full;final Locus locus;
	}

	private final String family;private final int flank;private final boolean panel;private int runCount;
	private final ObjectMap<String,Locus> physical=new ObjectMap<String,Locus>(String.class, Locus.class);
	private final ObjectMap<String,Group> byContig=new ObjectMap<String,Group>(String.class, Group.class);
	private final ArrayList<Locus> loci=new ArrayList<Locus>();private final ArrayList<Member> members=new ArrayList<Member>();
	public static final String CALL_HEADER="call_id\tseqid\tfamily\tstart1\tstop1\tstrand\tsource_outcome\tfull_outcome\texact_duplicate\tindependent_review";
	public static final String PANEL_CALL_HEADER="view\tscope\tcall_index\tgenome\tseqid\temitted_family\tstart1\tstop1\tstrand\tmodel\toutcome\ttruth_id\td5\td3";
	public static final String LOCUS_HEADER="review_id\tgenome\tassembly\tfasta\tseqid\tfamily\tstrand\tstart1\tstop1\tcontig_length\tlflank\trflank\tcontext_length\tstatus";
}
