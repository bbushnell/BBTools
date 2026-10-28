package prok;

import java.util.ArrayList;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.ReadWriteLock;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import idaligner.QuantumAligner;
import idaligner.ScrabbleAligner;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;
import structures.IntList;
import template.Accumulator;
import template.ThreadWaiter;

/** Exhaustive initial-verifier evidence on frozen nonedge caller windows.
 * Mapping curves reuse native index counts; they are not final caller P/R.
 * @author Brian, Raiden
 */
public final class NcrnaMappingOracle implements Accumulator<NcrnaMappingOracle.Worker> {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaMappingOracle.class, false);
		String input=null, windows=null, truth=null, ref=null, out=null, saved=null;int expected=-1, threads=1, minimum=109, quantum=120;float cut=.65f;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");
			final String key=arg.substring(0, eq).toLowerCase(Locale.ROOT), value=arg.substring(eq+1);
			if(key.equals("in")){input=value;}else if(key.equals("windows")){windows=value;}else if(key.equals("truth")){truth=value;}
			else if(key.equals("ref")){ref=value;}else if(key.equals("outprefix")){out=value;}
			else if(key.equals("scores")){saved=value;}
			else if(key.equals("expected")){expected=Integer.parseInt(value);}else if(key.equals("t")){threads=Integer.parseInt(value);}
			else if(key.equals("minlen")){minimum=Integer.parseInt(value);}else if(key.equals("quantum")){quantum=Integer.parseInt(value);}
			else if(key.equals("idpass")){cut=Float.parseFloat(value);}else{throw new IllegalArgumentException("Unknown option: "+key);}
		}
		require(input!=null && windows!=null && truth!=null && ref!=null && out!=null && expected>0 && threads>0
			&& minimum>0 && quantum>0 && cut>=0 && cut<=1, "Bind windows, CM-supported truth, library and the initial verifier settings");
		Shared.setThreads(threads);Shared.TRIM_READ_DESCRIPTION=false;
		final NcrnaMappingOracle tool=new NcrnaMappingOracle(ref, minimum, quantum, cut);
		tool.loadTruth(truth);tool.loadWindows(windows, input, expected);
		if(saved!=null){tool.loadScores(saved);}
		tool.process(out, threads);Shared.closeStream(pp.outstream);
	}
	private NcrnaMappingOracle(String ref, int minimum_, int quantum_, float cut_){
		library=TrnaConsensusBuilder.loadLibrary(ref);minimum=minimum_;quantum=quantum_;cut=cut_;
		require(library!=null && library.length>1, "This oracle measures the multi-model keep-best shortlist policy");
	}

	/*--------------------------------------------------------------*/
	/*----------------            Inputs            ----------------*/
	/*--------------------------------------------------------------*/

	private void loadTruth(final String path){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');final ObjectSet<String> unique=new ObjectSet<String>(String.class);
		try{
			RrnaResourceIO.header(in, NcrnaCalibrationPrep.TRUTH_HEADER, path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==5 && p.termEquals("cobalamin", 1), "Expected cobalamin full-model coordinate view");
				final Target t=new Target(p.parseString(0), p.parseInt(2)-1, p.parseInt(3)-1, p.termEquals("+", 4)?0:p.termEquals("-", 4)?1:-1);
				require(t.start>=0 && t.stop>=t.start && t.strand>=0 && unique.add(t.contig+'\t'+t.start+'\t'+t.stop+'\t'+t.strand),
					"Truth coordinates and physical identities must be unambiguous");targets.add(t);
			}
		}finally{require(!in.close(), "Truth read failed");}
		require(!targets.isEmpty(), "Oracle requires a nonempty labeled target view");
	}
	private void loadWindows(final String path, final String fasta, final int expected){
		final ObjectMap<String,Window> ids=new ObjectMap<String,Window>(String.class, Window.class);
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, NcrnaCalibrationPrep.WINDOW_HEADER, path);
			for(byte[] row; (row=in.nextLine())!=null; ){
				p.set(row);require(p.terms()==13, "Frozen window metadata schema changed");
				final Window w=new Window(p.parseString(0), p.parseString(2), p.parseInt(3), p.parseInt(4), p.parseInt(6), p.parseInt(7),
					p.parseString(8), p.parseInt(9), p.parseInt(10), p.parseFloat(11), p.termEquals("true", 12));
				require(!w.id.isEmpty() && !ids.contains(w.id) && w.start>=0 && w.stop>=w.start && w.stop<w.contigLength && (w.strand==0 || w.strand==1)
					&& (p.termEquals("true", 12) || p.termEquals("false", 12)) && w.edge==(w.start==0 || w.stop==w.contigLength-1)
					&& w.shortlisted>=-1 && w.shortlisted<=library.length && w.baselineModel>=-1 && w.baselineModel<library.length,
					"Frozen window geometry, physical-edge label and model indices must agree");
				for(int i=0; i<targets.size(); i++){
					final Target t=targets.get(i);if(!w.contig.equals(t.contig) || w.strand!=t.strand){continue;}
					require(t.stop<w.contigLength, "Target exceeds frozen contig length");
					final int a=w.strand==0?t.start:w.contigLength-1-t.stop, z=w.strand==0?t.stop:w.contigLength-1-t.start;
					if(a<=w.stop && z>=w.start){w.targets.add(i);}
				}
				ids.put(w.id, w);windows.add(w);
			}
		}finally{require(!in.close(), "Window metadata read failed");}
		require(windows.size()==expected, "Frozen window event denominator changed");
		RrnaResourceIO.stream(fasta, r->{
			final Window w=ids.get(RrnaResourceIO.first(r.id));
			require(w!=null && w.seq==null && r.bases!=null && r.length()==w.stop-w.start+1, "Each window needs its one exact extracted sequence: "+r.id);
			w.seq=r.bases;
		});
		for(Window w:windows){require(w.seq!=null, "Missing window sequence: "+w.id);}
	}
	/** Reuses complete ordered scores; native shared counts and baseline replay are still verified. */
	private void loadScores(final String path){
		final ByteFile in=RrnaResourceIO.open(path);final LineParser1 p=new LineParser1('\t');
		try{
			RrnaResourceIO.header(in, SCORE_HEADER, path);
			for(Window w:windows){
				if(w.edge || w.shortlisted<0){continue;}
				w.identity=new float[library.length];w.starts=new int[library.length];w.stops=new int[library.length];w.shared=new int[KS.length][library.length];
				for(int m=0; m<library.length; m++){
					final byte[] row=in.nextLine();require(row!=null, "Saved scores are incomplete for "+w.id);p.set(row);
					require(p.terms()==9 && p.termEquals(w.id, 0) && p.parseInt(1)==m, "Saved scores must retain each ordered window/model pair exactly once");
					w.identity[m]=p.parseFloat(2);w.starts[m]=p.parseInt(3);w.stops[m]=p.parseInt(4);
					for(int k=0; k<KS.length; k++){w.shared[k][m]=p.parseInt(5+k);require(w.shared[k][m]>=0, "Saved shared counts cannot be negative");}
				}
			}
			require(in.nextLine()==null, "Saved scores include extra or excluded window/model pairs");
		}finally{require(!in.close(), "Saved model score read failed");}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Thread Management      ----------------*/
	/*--------------------------------------------------------------*/

	private void process(final String out, final int threads){
		RrnaResourceIO.fresh(out, ".models.tsv.gz", ".windows.tsv", ".curves.tsv", ".losses.tsv", ".loss_summary.tsv");
		covered=new long[CONFIGS][(targets.size()+63)/64];oracleCovered=new long[(targets.size()+63)/64];
		recovery=new byte[targets.size()];
		final ArrayList<Worker> workers=new ArrayList<Worker>();
		for(int i=0; i<Math.min(threads, windows.size()); i++){workers.add(new Worker(this));}
		if(!ThreadWaiter.startAndWait(workers, this)){throw new IllegalStateException("Oracle worker failed; no reports written", failure);}
		require(calls+reused==tested*library.length, "Every considered window must have one new or reused score for every model");
		write(out);
	}
	@Override
	public void accumulate(final Worker w){
		assert(!w.isAlive()):"ThreadWaiter joins before publishing oracle results";
		if(!w.success && failure==null){failure=w.failure;}
		failed|=!w.success;calls+=w.calls;reused+=w.reused;tested+=w.tested;oracleEvents+=w.oracleEvents;
		for(int c=0; c<CONFIGS; c++){
			admitted[c]+=w.admitted[c];retainedEvents[c]+=w.retainedEvents[c];
			for(int i=0; i<covered[c].length; i++){covered[c][i]|=w.covered[c][i];}
		}
		for(int i=0; i<oracleCovered.length; i++){oracleCovered[i]|=w.oracleCovered[i];}
		for(int i=0; i<recovery.length; i++){recovery[i]=(byte)Math.max(recovery[i], w.recovery[i]);}
	}
	@Override
	public boolean success(){return !failed;}
	@Override
	public ReadWriteLock rwlock(){return null;}

	/** Counts a prefix of native count-desc/id-asc ordering, preserving keep-best. */
	public static int retainedCount(final int[] order, final int[] shared, final int top, final int minimum){
		assert(top>0 && minimum>=0):"The mapping grid must use valid native shortlist settings";
		int n=0;while(n<order.length && n<top && shared[order[n]]>=minimum){n++;}
		return n==0 && order.length>0?1:n;
	}
	private boolean matches(final Window w, final int model, final int target){
		if(!eligible(w, model)){return false;}
		final Target t=targets.get(target);
		final int a=w.strand==0?t.start:w.contigLength-1-t.stop, z=w.strand==0?t.stop:w.contigLength-1-t.start;
		final int start=w.start+w.starts[model], stop=w.start+w.stops[model];
		final int overlap=Math.min(stop, z)-Math.max(start, a)+1;
		return overlap>0 && 2L*overlap>=stop-start+1 && 2L*overlap>=z-a+1;
	}
	private boolean eligible(final Window w, final int m){return w.identity[m]>=cut && w.stops[m]-w.starts[m]+1>=minimum;}
	private static void mark(final long[] bits, final int index){bits[index>>>6]|=1L<<(index&63);}
	private static int count(final long[] bits){int n=0;for(long value:bits){n+=Long.bitCount(value);}return n;}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	static final class Worker extends Thread {
		Worker(final NcrnaMappingOracle owner_){
			owner=owner_;covered=new long[CONFIGS][owner.oracleCovered.length];oracleCovered=new long[owner.oracleCovered.length];
			recovery=new byte[owner.targets.size()];
			for(int k=0; k<KS.length; k++){indexes[k]=new TrnaKmerIndex(owner.library, KS[k], false, 0, 0, 0, 1);}
		}
		@Override
		public void run(){
			try{
				for(int first; (first=owner.next.getAndAdd(8))<owner.windows.size(); ){
					for(int i=first, end=Math.min(first+8, owner.windows.size()); i<end; i++){process(owner.windows.get(i));}
				}
				success=true;
			}catch(Throwable t){failure=t;owner.next.set(owner.windows.size());}
		}
		private void process(final Window w){
			if(w.edge || w.shortlisted<0){w.status=w.edge?"PHYSICAL_EDGE_SEPARATE":"PREINDEX_GATE";return;}
			w.status="SCORED";tested++;
			final int models=owner.library.length;final boolean saved=w.identity!=null;
			if(!saved){w.identity=new float[models];w.starts=new int[models];w.stops=new int[models];w.shared=new int[KS.length][models];}
			final int[] pos=new int[4];boolean any=false;
			for(int m=0; m<models; m++){
				if(saved){reused++;}
				else{
					if(w.seq.length>owner.quantum){w.identity[m]=QuantumAligner.alignStatic(owner.library[m], w.seq, pos);w.starts[m]=pos[0];w.stops[m]=pos[1];}
					else{w.identity[m]=ScrabbleAligner.alignStatic(w.seq, owner.library[m], null);w.starts[m]=0;w.stops[m]=w.seq.length-1;}calls++;
				}
				require(Float.isFinite(w.identity[m]) && w.identity[m]>=0 && w.identity[m]<=1, "Nonfractional oracle identity");
				if(w.identity[m]>=owner.cut){require(w.starts[m]>=0 && w.stops[m]>=w.starts[m] && w.stops[m]<w.seq.length,
					"Passing initial alignments must report a real span inside their frozen query window");}
				if(owner.eligible(w, m)){
					any=true;
					for(int j=0; j<w.targets.size; j++){final int t=w.targets.get(j);if(owner.matches(w, m, t)){mark(oracleCovered, t);}}
				}
			}
			if(any){oracleEvents++;}
			for(int k=0; k<KS.length; k++){
				final int[] order=indexes[k].shortlist(w.seq, models);
				for(int m=0; m<models; m++){
					final int count=indexes[k].lastSharedCount(m);
					if(saved){require(w.shared[k][m]==count, "Saved counts must match the actual native index on frozen sequences");}else{w.shared[k][m]=count;}
				}
				if(KS[k]==7){checkBaseline(w, order);classifyLosses(w, order);}
				for(int n=0; n<TOPS.length; n++){
					for(int h=0; h<HITS.length; h++){
						final int c=(k*TOPS.length+n)*HITS.length+h, size=retainedCount(order, w.shared[k], TOPS[n]==0?models:TOPS[n], HITS[h]);
						admitted[c]+=size;boolean pass=false;
						for(int j=0; j<size; j++){
							final int m=order[j];if(!owner.eligible(w, m)){continue;}pass=true;
							for(int z=0; z<w.targets.size; z++){final int t=w.targets.get(z);if(owner.matches(w, m, t)){mark(covered[c], t);}}
						}
						if(pass){retainedEvents[c]++;}
					}
				}
			}
		}
		private void classifyLosses(final Window w, final int[] order){
			assert(w.identity!=null):"Cap attribution requires exhaustive initial-verifier scores";
			final int[] rank=new int[owner.library.length];for(int i=0; i<order.length; i++){rank[order[i]]=i+1;}
			for(int j=0; j<w.targets.size; j++){
				final int t=w.targets.get(j);byte best=1;
				for(int m=0; m<rank.length; m++){
					if(!owner.matches(w, m, t)){continue;}
					final byte state=(byte)(rank[m]==0?2:rank[m]>30?3:4);best=(byte)Math.max(best, state);
				}
				recovery[t]=(byte)Math.max(recovery[t], best);
			}
		}
		private void checkBaseline(final Window w, final int[] order){
			final int size=retainedCount(order, w.shared[1], 30, 1);
			require(size==w.shortlisted, "Rebuilt baseline shortlist size differs for "+w.id);
			if(!(w.outcome.equals("ACCEPT_IDPASS") || w.outcome.equals("REJECT_ID") || w.outcome.equals("REJECT_LEN"))){return;}
			float best=0;int model=-1;
			for(int i=0; i<size; i++){final int m=order[i];if(w.identity[m]>best){best=w.identity[m];model=m;}}
			require(Float.floatToIntBits(best)==Float.floatToIntBits(w.baselineId) && model==w.baselineModel,
				"Oracle must replay the actual initial verifier before deriving new mapping curves: "+w.id);
		}
		final NcrnaMappingOracle owner;final TrnaKmerIndex[] indexes=new TrnaKmerIndex[KS.length];
		final long[][] covered;final long[] oracleCovered, admitted=new long[CONFIGS], retainedEvents=new long[CONFIGS];
		final byte[] recovery;
		long calls, reused, tested, oracleEvents;boolean success;Throwable failure;
	}
	private static final class Target {
		Target(String c, int a, int z, int s){contig=c;start=a;stop=z;strand=s;}
		final String contig;final int start, stop, strand;
	}
	private static final class Window {
		Window(String id_, String c, int s, int length, int a, int z, String outcome_, int shortlisted_, int model, float identity_, boolean edge_){
			id=id_;contig=c;strand=s;contigLength=length;start=a;stop=z;outcome=outcome_;shortlisted=shortlisted_;baselineModel=model;baselineId=identity_;edge=edge_;
		}
		final String id, contig, outcome;final int strand, contigLength, start, stop, shortlisted, baselineModel;final float baselineId;final boolean edge;
		final IntList targets=new IntList();byte[] seq;String status;float[] identity;int[] starts, stops;int[][] shared;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Outputs           ----------------*/
	/*--------------------------------------------------------------*/

	private void write(final String prefix){
		final ByteStreamWriter models=RrnaResourceIO.writer(prefix+".models.tsv.gz"), events=RrnaResourceIO.writer(prefix+".windows.tsv"), curves=RrnaResourceIO.writer(prefix+".curves.tsv");
		final ByteBuilder b=new ByteBuilder();final int oracleLoci=count(oracleCovered);
		try{
			models.println(SCORE_HEADER);
			events.println("window\tstatus\ttruth_overlapping\tbaseline_outcome");
			for(Window w:windows){
				events.print(b.clear().append(w.id).tab().append(w.status).tab().append(w.targets.size).tab().append(w.outcome).nl());
				if(w.identity==null){continue;}
				for(int m=0; m<library.length; m++){
					b.clear().append(w.id).tab().append(m).tab().appendSlow(w.identity[m]).tab().append(w.starts[m]).tab().append(w.stops[m]);
					for(int k=0; k<KS.length; k++){b.tab().append(w.shared[k][m]);}models.print(b.nl());
				}
			}
			curves.println("index_k\ttop_n\tmin_hits\tall_window_events\ttested_nonedge_events\toracle_passing_events\tretained_passing_events\tfull_truth_loci\toracle_loci\tretained_oracle_loci\tshortlisted_models\tmean_shortlist\tnew_alignments\tinterpretation");
			for(int k=0; k<KS.length; k++){
				for(int n=0; n<TOPS.length; n++){
					for(int h=0; h<HITS.length; h++){
						final int c=(k*TOPS.length+n)*HITS.length+h;
						curves.print(b.clear().append(KS[k]).tab().append(TOPS[n]==0?library.length:TOPS[n]).tab().append(HITS[h]).tab().append(windows.size()).tab().append(tested)
							.tab().append(oracleEvents).tab().append(retainedEvents[c]).tab().append(targets.size()).tab().append(oracleLoci).tab().append(count(covered[c]))
							.tab().append(admitted[c]).tab().append(tested==0?0:admitted[c]/(double)tested, 6).tab().append(calls)
							.append("\tINITIAL_VERIFIER_POTENTIAL_NOT_CALLER_RECALL\n"));
					}
				}
			}
		}finally{RrnaResourceIO.finish(models, events, curves);}
		final ByteStreamWriter losses=RrnaResourceIO.writer(prefix+".losses.tsv"), lossSummary=RrnaResourceIO.writer(prefix+".loss_summary.tsv");
		int supported=0, retained=0;final int[] categories=new int[LOSS_STATUS.length];
		try{
			losses.println("seqid\tstart1\tstop1\tstrand\tk7_min1_top30_status");
			for(int i=0; i<targets.size(); i++){
				final Target t=targets.get(i);supported+=recovery[i]>=2?1:0;retained+=recovery[i]==4?1:0;
				categories[recovery[i]]++;
				losses.print(b.clear().append(t.contig).tab().append(t.start+1).tab().append(t.stop+1).tab().append(t.strand==0?'+':'-').tab().append(LOSS_STATUS[recovery[i]]).nl());
			}
			lossSummary.println("status\tloci");
			for(int i=0; i<categories.length; i++){lossSummary.print(b.clear().append(LOSS_STATUS[i]).tab().append(categories[i]).nl());}
		}finally{RrnaResourceIO.finish(losses, lossSummary);}
		require(supported==oracleLoci && retained==count(covered[(TOPS.length+3)*HITS.length]), "Per-locus loss partition must reproduce oracle and baseline coverage");
		System.err.println("NCRNA_MAPPING_ORACLE_PASS events="+windows.size()+" tested="+tested+" models="+library.length+" alignments="+calls+" reused="+reused+" oracle_loci="+oracleLoci);
	}
	private static void require(final boolean ok, final String why){if(!ok){throw new IllegalArgumentException(why);}}
	private final byte[][] library;private final int minimum, quantum;private final float cut;
	private final ArrayList<Window> windows=new ArrayList<Window>();private final ArrayList<Target> targets=new ArrayList<Target>();
	private final AtomicInteger next=new AtomicInteger();private boolean failed;private Throwable failure;
	private long calls, reused, tested, oracleEvents;private long[][] covered;private long[] oracleCovered;private byte[] recovery;
	private final long[] admitted=new long[CONFIGS], retainedEvents=new long[CONFIGS];
	private static final int[] KS={6,7,8,9}, TOPS={1,3,10,30,50,100,0}, HITS={1,3,5,10,20,40,80};
	private static final String SCORE_HEADER="window\tmodel\tidentity\taligned_start0\taligned_stop0\tshared6\tshared7\tshared8\tshared9";
	private static final String[] LOSS_STATUS={"NO_SCORED_OVERLAPPING_WINDOW", "NO_ACCEPTABLE_ALIGNED_MODEL", "ACCEPTABLE_MODEL_NOT_INDEXED", "CAP_EVICTED", "RETAINED_TOP30"};
	private static final int CONFIGS=KS.length*TOPS.length*HITS.length;
}
