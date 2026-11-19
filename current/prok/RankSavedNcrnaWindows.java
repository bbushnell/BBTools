package prok;

import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

import cardinality.DynamicDemiLog;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlSavedData.Edge;
import prok.NcrnaDdlSavedData.Parent;
import prok.NcrnaDdlSavedData.Selected;
import prok.NcrnaDdlSavedData.Window;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;
import template.ThreadWaiter;

/** Saved-window forward/plain/middle-masked grid with exact canonical controls.
 * Uses the ByteFile/Streamer templates' parser, native I/O and worker lifecycle.
 * The owner loads the small saved window cache before workers start; workers
 * claim coarse blocks, own all scoring/minimum scratch, and merge after joining.
 * No genome, GeneCaller, alignment or new truth-admission path is present.
 * @author Brian, Raiden
 */
public final class RankSavedNcrnaWindows{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, RankSavedNcrnaWindows.class, false);
		final int previous=DynamicDemiLog.exponentBits();
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			final Config c=new Config();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(a.equals("panels")){c.panels=b;}
				else if(a.equals("selection")){c.selection=b;}
				else if(a.equals("members")){c.members=b;}
				else if(a.equals("ledger")){c.ledger=b;}
				else if(a.equals("assignments")){c.assignments=b;}
				else if(a.equals("census")){c.census=b;}
				else if(a.equals("summary")){c.summary=b;}
				else if(a.equals("ref")){c.references=b;}
				else if(a.equals("expectedgroups")){c.groups=Parse.parseIntKMG(b);}
				else if(a.equals("expectedmodels")){c.models=Parse.parseIntKMG(b);}
				else if(a.equals("t") || a.equals("threads")){c.threads=Parse.parseIntKMG(b);}
				else{require(parser.parse(arg, a, b), "Unknown saved-grid parameter: "+arg);}
			}
			require(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use panels= selection= out=fresh_prefix; no genome input, overwrite, append or sampling");
			c.output=parser.out1;
			DynamicDemiLog.setExponent(4);
			run(c);
			pp.outstream.println("NCRNA_SAVED_FORWARD_MASKED_GRID_COMPLETE_NOT_CALLER_ADOPTED");
		}finally{
			DynamicDemiLog.setExponent(previous);
			Shared.closeStream(pp.outstream);
		}
	}

	public static void run(final Config c){
		require(c!=null && c.panels!=null && c.selection!=null && c.output!=null && c.members!=null
			&& c.ledger!=null && c.assignments!=null && c.census!=null && c.summary!=null
			&& c.references!=null && c.models>0 && c.groups>0 && c.threads>0,
			"Declare the complete saved selection, frozen member construction and positive worker/model/group counts");
		RrnaResourceIO.fresh(c.output, ".methods.tsv", ".summary.tsv", ".loci.tsv.gz", ".orientation.tsv", ".last_losses.tsv", ".work.tsv");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions && DynamicDemiLog.exponentBits()==4 && !DynamicDemiLog.blacklistExists(),
			"Saved ranking requires assertions, native exponent4 and no process-global DDL blacklist");
		Shared.TRIM_READ_DESCRIPTION=false;
		final RankSavedNcrnaWindows tool=new RankSavedNcrnaWindows(c);
		tool.process();
	}

	private RankSavedNcrnaWindows(final Config c_){c=c_;}

	/*--------------------------------------------------------------*/
	/*----------------         Owner Prepasses      ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		final long start=cpu();
		data=NcrnaDdlSavedData.load(c.panels, c.selection, c.models, c.groups);
		windows=data.windows.toArray(new Window[data.windows.size()]);
		geometries=geometries();
		final NcrnaDdlMembers members=members();
		final long profileStart=cpu();
		profiles=new NcrnaDdlProfile[geometries.length][];
		for(int g=0; g<geometries.length; g++){
			final Geometry geometry=geometries[g];
			profiles[g]=NcrnaDdlModelProfiles.build(members, geometry.k, geometry.width, geometry.canonical, geometry.middle);
			profileKeys+=(long)geometry.width*c.models*Character.BYTES;
			profileSupport+=(long)geometry.width*c.models*Integer.BYTES;
		}
		profileCpu=cpu()-profileStart;
		final int workers=Math.min(c.threads, Math.max(1, windows.length));
		final ArrayList<ProcessThread> threads=new ArrayList<ProcessThread>(workers);
		for(int i=0; i<workers; i++){threads.add(new ProcessThread(this, i));}
		next.set(Math.min(workers, windows.length));
		ThreadWaiter.startAndWait(threads);
		if(failure.get()!=null){throw new IllegalStateException("Saved ranking failed; no result tables were written", failure.get());}
		minima=new int[data.parents.size()*METHODS];
		Arrays.fill(minima, -1);
		for(final ProcessThread worker:threads){
			require(worker.success, "A saved-grid worker must finish normally before its minima are used");
			for(int i=0; i<minima.length; i++){minima[i]=NcrnaDdlSavedData.minimum(minima[i], worker.minima[i]);}
			workerCpu+=worker.cpuNanos;
			scoredWindows+=worker.scored;
			parityTuples+=worker.parity;
			if(worker.scored>0){scoringWorkers++;}
		}
		for(final Parent parent:data.parents){
			for(int control=0; control<NcrnaDdlSavedData.CONTROLS; control++){
				final int method=FORWARD_METHODS+control;
				require(minima[parent.index*METHODS+method]==parent.controlMin[control],
					"Recomputed canonical locus minima must exactly match the accepted saved controls");
			}
		}
		report();
		try(Output work=new Output(c.output+".work.tsv")){
			final ByteBuilder b=new ByteBuilder();
			work.writer.print(b.append("item\tvalue\nselected_inputs\t").append(data.selected.size())
				.append("\nselected_targets\t").append(data.parents.size()).append("\nselected_windows\t").append(windows.length)
				.append("\nscored_windows\t").append(scoredWindows).append("\ncanonical_event_parent_score_tuples_checked\t").append(parityTuples)
				.append("\nworkers_started\t").append(workers).append("\nworkers_with_scored_windows\t").append(scoringWorkers)
				.append("\nprofile_key_payload_bytes\t").append(profileKeys).append("\nprofile_support_payload_bytes\t").append(profileSupport)
				.append("\nprofile_heap_requirement\tUNMEASURED_NOT_PAYLOAD_BYTES\nprofile_owner_cpu_ns\t").append(profileCpu)
				.append("\nworker_cpu_ns\t").append(workerCpu).append("\nowner_cpu_ns\t").append(cpu()-start)
				.append("\ncpu_scope\tOWNER_AND_WORKERS_ONLY_EXCLUDES_NATIVE_IO_AND_GC\nresult\tSAVED_RANK_MINIMA_NOT_CALLER_ALIGNMENTS_OR_ADOPTION\n"));
		}
	}

	private NcrnaDdlMembers members(){
		final ArrayList<String> names=new ArrayList<String>();
		final ArrayList<byte[]> bases=new ArrayList<byte[]>();
		NcrnaDdlSavedData.sequences(c.references, read->{
			names.add(RrnaResourceIO.first(read.id));
			bases.add(read.bases);
		});
		require(names.size()==c.models, "Frozen reference order must retain every declared original model");
		return NcrnaDdlMembers.load(c.members, c.ledger, c.assignments, c.census, c.summary,
			names.toArray(new String[names.size()]), bases.toArray(new byte[bases.size()][]));
	}

	/*--------------------------------------------------------------*/
	/*----------------          Grid Reports        ----------------*/
	/*--------------------------------------------------------------*/

	private void report(){
		final long[][][] hist=new long[3][METHODS+1][c.models+1];
		final long[][] sums=new long[3][METHODS+1];
		final int[] common=new int[3];
		final boolean[] baselineReach=new boolean[data.selected.size()];
		final boolean[][] candidateReach=new boolean[METHODS][data.selected.size()];
		try(Output methods=new Output(c.output+".methods.tsv"); Output loci=new Output(c.output+".loci.tsv.gz");
			Output stats=new Output(c.output+".summary.tsv"); Output orientation=new Output(c.output+".orientation.tsv");
			Output losses=new Output(c.output+".last_losses.tsv")){
			final ByteBuilder b=new ByteBuilder();
			methods.writer.println("method\tk\tbuckets\torientation\tmiddle_mask\tscore\trole");
			methods.writer.println("0\t9\t0\tFORWARD_NATIVE_INDEX\tfalse\tPOSITIVE_SHARED_KEYS\tSAVED_BASELINE");
			for(int g=0; g<geometries.length; g++){
				final Geometry geometry=geometries[g];
				for(int score=0; score<2; score++){
					methods.writer.print(b.clear().append(1+2*g+score).tab().append(geometry.k).tab().append(geometry.width)
						.tab().append(geometry.canonical ? "CANONICAL" : "FORWARD").tab().append(geometry.middle)
						.tab().append(SCORE_NAMES[score]).tab().append(geometry.canonical ? "EXACT_SAVED_CONTROL" : "NEW_GRID").nl());
				}
			}
			loci.writer.println("panel\tinput_slot\tinput_id\tparent\tseen_stratum\tmethod\tbest_rank\tsource_disposition\tinterpretation");
			for(final Parent parent:data.parents){
				final int offset=parent.index*METHODS;
				if(parent.nativeRank>0){common[parent.stratum]++;}
				if(parent.nativeRank>0 && parent.nativeRank<=50){baselineReach[parent.input.index]=true;}
				for(int method=-1; method<METHODS; method++){
					final int rank=method<0 ? parent.nativeRank : minima[offset+method];
					if(method>=0 && rank>0 && rank<=50){candidateReach[method][parent.input.index]=true;}
					if(parent.nativeRank>0){
						require(rank>0 && rank<=c.models, "Every common locus must retain a positive rank in each new/control method");
						hist[parent.stratum][method+1][rank]++;
						sums[parent.stratum][method+1]+=rank;
					}
					loci.writer.print(b.clear().append(parent.input.panel.nameBytes).tab().append(parent.input.slot)
						.tab().append(parent.input.idBytes).tab().append(parent.idBytes).tab().append(NcrnaDdlSavedData.STRATA[parent.stratum])
						.tab().append(method+1).tab().append(rank).tab().append(parent.reason)
						.append("\tSAVED_LOCUS_MINIMUM_NOT_CALLER_WORK\n"));
				}
			}
			stats.writer.println("seen_stratum\tmethod\tcommon_loci\trank_sum\tmean_rank\twithin10\twithin50\twithin100\tp90\tp95\tp99\tinterpretation");
			for(int stratum=0; stratum<3; stratum++){
				for(int method=0; method<=METHODS; method++){
					final long[] h=hist[stratum][method];
					b.clear().append(NcrnaDdlSavedData.STRATA[stratum]).tab().append(method).tab().append(common[stratum]).tab().append(sums[stratum][method]).tab();
					if(common[stratum]==0){b.append("NA");}else{b.append((double)sums[stratum][method]/common[stratum], 9);}
					stats.writer.print(b.tab().append(covered(h, 10)).tab().append(covered(h, 50)).tab().append(covered(h, 100))
						.tab().append(RankNcrnaWindows.quantile(h, .90)).tab().append(RankNcrnaWindows.quantile(h, .95))
						.tab().append(RankNcrnaWindows.quantile(h, .99)).append("\tSAVED_LOCUS_MINIMUM_NOT_CALLER_ALIGNMENTS\n"));
				}
			}
			orientation(orientation, sums, common);
			losses.writer.println("method\tpanel\tinput_slot\tinput_id\tfasta\tinterpretation");
			for(int method=0; method<METHODS; method++){
				for(final Selected input:data.selected){
					if(baselineReach[input.index] && !candidateReach[method][input.index]){
						losses.writer.print(b.clear().append(method+1).tab().append(input.panel.nameBytes).tab().append(input.slot)
							.tab().append(input.idBytes).tab().append(input.pathBytes)
							.append("\tLAST_SAVED_RANK_REACHABLE_LOCUS_LOST_AT50_NOT_CALLER_DETECTION\n"));
					}
				}
			}
		}
	}

	private void orientation(final Output out, final long[][] sums, final int[] common){
		out.writer.println("seen_stratum\tk\tbuckets\tscore\tcommon_loci\tforward_rank_sum\tcanonical_rank_sum\tforward_better\ttied\tforward_worse");
		final ByteBuilder b=new ByteBuilder();
		for(int control=0; control<NcrnaDdlSavedData.CONTROLS; control++){
			final Geometry canonical=geometries[FORWARD_GEOMETRIES+control/2];
			final int widthIndex=(control/2)%4;
			final int forward=2*((canonical.k-7)*8+widthIndex)+(control&1);
			final int comparator=FORWARD_METHODS+control;
			for(int stratum=0; stratum<3; stratum++){
				int better=0, tied=0, worse=0;
				for(final Parent parent:data.parents){
					if(parent.stratum!=stratum || parent.nativeRank<0){continue;}
					final int a=minima[parent.index*METHODS+forward], c=minima[parent.index*METHODS+comparator];
					if(a<c){better++;}else if(a==c){tied++;}else{worse++;}
				}
				out.writer.print(b.clear().append(NcrnaDdlSavedData.STRATA[stratum]).tab().append(canonical.k).tab().append(canonical.width)
					.tab().append(SCORE_NAMES[control&1]).tab().append(common[stratum]).tab().append(sums[stratum][forward+1])
					.tab().append(sums[stratum][comparator+1]).tab().append(better).tab().append(tied).tab().append(worse).nl());
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Worker Ownership     ----------------*/
	/*--------------------------------------------------------------*/

	static final class ProcessThread extends Thread{
		ProcessThread(final RankSavedNcrnaWindows owner_, final int index_){
			owner=owner_;
			index=index_;
			minima=new int[owner.data.parents.size()*METHODS];
			Arrays.fill(minima, -1);
			scoring=new NcrnaDdlRanks[owner.geometries.length];
			for(int g=0; g<scoring.length; g++){
				final Geometry geometry=owner.geometries[g];
				scoring[g]=new NcrnaDdlRanks(owner.profiles[g], geometry.k, geometry.canonical, geometry.middle);
			}
			acceptable=new boolean[owner.c.models];
		}
		@Override
		public void run(){
			final long start=cpu();
			try{
				// One initial window each; remaining work uses coarse atomic claims.
				// This also guarantees genuine cross-worker fixture coverage.
				if(index<owner.windows.length){score(owner.windows[index]);}
				for(int first=owner.next.getAndAdd(8); first<owner.windows.length; first=owner.next.getAndAdd(8)){
					if(owner.failure.get()!=null){return;}
					for(int i=first; i<Math.min(first+8, owner.windows.length); i++){score(owner.windows[i]);}
				}
				success=true;
			}catch(RuntimeException | Error invalid){owner.failure.compareAndSet(null, invalid);}
			finally{cpuNanos=cpu()-start;}
		}
		private void score(final Window window){
			boolean useful=false;
			for(final Edge edge:window.edges){useful|=edge.unique>0;}
			if(useful){
				for(final NcrnaDdlRanks ranks:scoring){ranks.score(window.bases);}
				scored++;
			}
			for(final Edge edge:window.edges){
				Arrays.fill(acceptable, false);
				for(int model=0; model<acceptable.length; model++){
					acceptable[model]=(edge.uniqueModels[model>>>6]&(1L<<(model&63)))!=0;
				}
				for(int g=0; g<scoring.length; g++){
					for(int score=0; score<2; score++){
						final int method=2*g+score;
						if(edge.unique>0){scoring[g].best(SCORES[score], acceptable, tuple);}
						else{Arrays.fill(tuple, -1);}
						if(g>=FORWARD_GEOMETRIES){
							final int offset=(method-FORWARD_METHODS)*4;
							for(int i=0; i<4; i++){
								require(tuple[i]==edge.controls[offset+i], "Canonical saved event rank/tie/model parity failed: panel="
									+window.input.panel.name+", event="+window.key+", control="+(method-FORWARD_METHODS)+", field="+i);
							}
							parity++;
						}
						final int index=edge.parent.index*METHODS+method;
						minima[index]=NcrnaDdlSavedData.minimum(minima[index], tuple[0]);
					}
				}
			}
		}
		private final RankSavedNcrnaWindows owner;
		private final int index;
		private final NcrnaDdlRanks[] scoring;
		private final boolean[] acceptable;
		private final int[] tuple=new int[4];
		final int[] minima;
		boolean success;
		long cpuNanos, scored, parity;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Configuration       ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Config{
		public String panels, selection, output, members, ledger, assignments, census, summary, references;
		public int models=508, groups=3153, threads=1;
	}
	private static final class Geometry{
		Geometry(final int k_, final int width_, final boolean canonical_, final boolean middle_){
			assert(k_>=7 && k_<=11 && Integer.bitCount(width_)==1 && !(canonical_ && middle_))
				: "The reviewed grid binds member/query geometry identically, with canonical plain controls only";
			k=k_;
			width=width_;
			canonical=canonical_;
			middle=middle_;
		}
		final int k, width;
		final boolean canonical, middle;
	}
	private static Geometry[] geometries(){
		final Geometry[] result=new Geometry[GEOMETRIES];
		int next=0;
		for(int k=7; k<=11; k++){
			for(final boolean middle:new boolean[]{false, true}){
				for(final int width:WIDTHS){result[next++]=new Geometry(k, width, false, middle);}
			}
		}
		for(final int k:new int[]{9, 11}){
			for(final int width:WIDTHS){result[next++]=new Geometry(k, width, true, false);}
		}
		assert(next==result.length) : "Forty forward geometries and eight canonical geometries define the fixed saved grid";
		return result;
	}
	private static long covered(final long[] histogram, final int limit){
		long count=0;
		for(int i=1; i<=limit && i<histogram.length; i++){count+=histogram[i];}
		return count;
	}
	private static long cpu(){
		final ThreadMXBean bean=ManagementFactory.getThreadMXBean();
		require(bean.isCurrentThreadCpuTimeSupported(), "Saved-grid diagnostic accounting requires a current-thread CPU clock");
		if(!bean.isThreadCpuTimeEnabled()){bean.setThreadCpuTimeEnabled(true);}
		return bean.getCurrentThreadCpuTime();
	}
	private static void require(final boolean valid, final String why){NcrnaDdlSavedData.require(valid, why);}

	private final Config c;
	private NcrnaDdlSavedData data;
	private Window[] windows;
	private Geometry[] geometries;
	private NcrnaDdlProfile[][] profiles;
	private int[] minima;
	private int scoringWorkers;
	private long profileKeys, profileSupport, profileCpu, workerCpu, scoredWindows, parityTuples;
	private final AtomicInteger next=new AtomicInteger();
	private final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
	private static final int FORWARD_GEOMETRIES=40, GEOMETRIES=48, FORWARD_METHODS=80, METHODS=96;
	private static final int[] WIDTHS={256,512,1024,2048}, SCORES={NcrnaDdlRanks.SSU,NcrnaDdlRanks.ALL_KEYS};
	private static final String[] SCORE_NAMES={"SSU","RAW_ALL_QUERY_KEYS"};
}
