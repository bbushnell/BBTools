package prot;

import java.io.PrintStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicReference;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;

import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.ByteStreamWriter;
import shared.KillSwitch;
import stream.Read;
import stream.Streamer;
import stream.StreamerFactory;
import structures.ListNum;
import structures.ByteBuilder;
import template.Accumulator;
import template.ThreadWaiter;
import prot.FamilyPanelRunner.Panel;
import prot.FamilyPanelRunner.PanelCounts;

/** One worker pool scores three ranking arms jointly on complete labeled panels.
 * Shared per-query scoring is supplied by a worker-local factory; no policy is selected here.
 * @author Sayu
 */
final class FamilyPairedPanelRunner implements Accumulator<FamilyPairedPanelRunner.ProcessThread>{

	FamilyPairedPanelRunner(FamilyTailCascade.LoadedClassifier classifier_,Panel[] panels_,int threads_,int inputThreads_,
			ScorerFactory factory_,PrintStream outstream_,String label_){
		require(classifier_!=null && classifier_.familyCount()>0 && panels_!=null && panels_.length>0,
			"Panel execution requires a loaded live classifier and nonempty panel roster");
		require(threads_>0 && inputThreads_>0 && outstream_!=null && label_!=null,
			"Panel execution needs positive worker counts and a progress destination");
		classifier=classifier_;panels=panels_;threads=threads_;inputThreads=inputThreads_;
		factory=factory_;outstream=outstream_;label=label_;
		require(factory!=null,"Paired evaluation requires a worker-local scorer factory");
		counts=new FamilyPairedCounts[panels.length];
		for(int i=0;i<counts.length;i++){counts[i]=new FamilyPairedCounts(classifier.familyCount());}
	}

	/** The caller pins inputs and publishes output; this method only scores and reduces. */
	void process(){
		require(!started,"A panel runner represents one evaluation; use a fresh runner for each round");started=true;
		final Streamer[] inputs=new Streamer[panels.length];
		for(int i=0;i<inputs.length;i++){
			final FileFormat format=FileFormat.testInput(panels[i].path.toString(),FileFormat.FASTA,null,true,true);
			require(format.fasta() || format.fastq(),"Protein panels must be FASTA or FASTQ: "+panels[i].path);
			inputs[i]=StreamerFactory.makeStreamer(format,null,true,-1,false,true,inputThreads);
		}
		final ArrayList<ProcessThread> workers=new ArrayList<ProcessThread>(threads);
		for(int i=0;i<threads;i++){workers.add(new ProcessThread(i,inputs));}
		try{
			for(Streamer input:inputs){input.start();}
			if(!ThreadWaiter.startAndWait(workers,this)){errorState=true;}
		}catch(Throwable problem){
			problem.printStackTrace(outstream);
			KillSwitch.kill("Family panel execution could not start or join its scoring workers");
		}
		for(Streamer input:inputs){errorState|=ReadWrite.closeStream(input);}
		if(failure.get()!=null){
			throw new IllegalStateException(
				"Family panel execution failed; no completed report is available",failure.get());
		}
		require(!errorState,"Family panel input or worker failed");
		long total=0;
		for(int i=0;i<counts.length;i++){
			counts[i].validate(panels[i].name);total=Math.addExact(total,counts[i].queries());
			require(panels[i].expectedQueries<0 || panels[i].expectedQueries==counts[i].queries(),
				"Panel "+panels[i].name+" observed count differs from its declared population");
		}
		require(total==progress.get(),"Worker totals disagree with processed input batches");
		complete=true;
	}

	/** Every arm uses the same observed panel denominator and explicit margin scale. */
	void writeSummary(Path path,double[] scales){
		require(complete && scales!=null && scales.length==3,"Completed paired results need the three frozen scales");
		final ByteBuilder out=new ByteBuilder("panel\tclass\tarm\tmargin_scale\tqueries\taccepted"+
			"\trejected\tTPR\tFNR\tFPR\tfalse_accepts_per_million");
		for(String stage:FamilyPanelRunner.METRIC_STAGES){out.tab().append(stage);}out.nl();
		for(int i=0;i<panels.length;i++){
			for(int a=0;a<3;a++){
				final Panel p=panels[i];final PanelCounts c=counts[i].arms[a];final double rate=c.accepted/(double)c.queries;
				out.append(p.name).tab().append(p.positive ? "positive" : "negative").tab().append(FamilyPairedCounts.ARM_NAMES[a])
					.tab().append(Double.toString(scales[a])).tab().append(c.queries).tab()
					.append(c.accepted).tab().append(c.queries-c.accepted).tab();
				if(p.positive){out.appendSlow(rate).tab().appendSlow(1-rate).append("\tNA\tNA");}
				else{out.append("NA\tNA\t").appendSlow(rate).tab().appendSlow(rate*1000000);}
				for(long n:c.stages){out.tab().append(n);}out.nl();
			}
		}
		write(path,out);
	}

	void writePairs(Path transitions,Path oracle){
		require(complete,"Paired outputs require all workers to finish");
		final ByteBuilder paired=new ByteBuilder(FamilyPairedCounts.TRANSITION_HEADER);
		final ByteBuilder masks=new ByteBuilder("panel\tqueries\tunion_alignment_pairs\toracle_membership_mask\tcount\n");
		for(int p=0;p<panels.length;p++){
			paired.append(counts[p].transitions(panels[p].name));
			for(int mask=0;mask<8;mask++){
				masks.append(panels[p].name).tab().append(counts[p].queries()).tab().append(counts[p].unionPairs())
					.tab().append(mask).tab().append(counts[p].oracleCount(mask)).nl();
			}
		}
		write(transitions,paired);write(oracle,masks);
	}

	/** Terminal-winner counts retain original physical ranks, including zero counts. */
	void writeFamilies(Path path){
		require(complete,"Per-family results require a completed paired evaluation");
		final ByteStreamWriter writer=new ByteStreamWriter(path.toString(),false,false,false);writer.start();
		try{
			writer.println("panel\tarm\troster_rank\tfamily_id\tnegative_queries\tfalse_accepts\tfalse_accepts_per_million");
			final ByteBuilder line=new ByteBuilder();
			for(int p=0;p<panels.length;p++){
				if(panels[p].positive){continue;}
				for(int a=0;a<3;a++){
					final PanelCounts c=counts[p].arms[a];
					for(int f=0;f<classifier.familyCount();f++){
						line.clear().append(panels[p].name).tab().append(FamilyPairedCounts.ARM_NAMES[a]).tab()
							.append(classifier.familyRankAt(f)).tab().append(classifier.familyIdAt(f)).tab().append(c.queries)
							.tab().append(c.winners[f]).tab().appendSlow(c.winners[f]*1000000.0/c.queries).nl();writer.print(line);
					}
				}
			}
		}finally{require(!writer.poisonAndWait(),"Paired family writer failed");}
	}

	@Override public void accumulate(ProcessThread worker){
		assert(worker.countsT.length==counts.length) : "Each worker must use the same ordered panel roster";
		for(int i=0;i<counts.length;i++){counts[i].add(worker.countsT[i]);}
		errorState|=!worker.success;
	}
	@Override public boolean success(){return !errorState && failure.get()==null;}
	@Override public ReadWriteLock rwlock(){return lock;}

	/** Streamers distribute batches; scratch and observations never cross worker threads. */
	final class ProcessThread extends Thread{
		ProcessThread(int id,Streamer[] inputs_){
			super("FamilyPanel-"+id);inputs=inputs_;
			countsT=new FamilyPairedCounts[panels.length];
			for(int i=0;i<countsT.length;i++){countsT[i]=new FamilyPairedCounts(classifier.familyCount());}
		}
		@Override public void run(){
			try{
				final Scorer scorer=factory.newScorer();
				for(int panel=0;panel<inputs.length;panel++){
					for(ListNum<Read> batch=inputs[panel].nextList();batch!=null;batch=inputs[panel].nextList()){
						try{
							if(failure.get()==null){
								currentCounts=countsT[panel];process(batch,scorer);
							}
						}catch(Throwable problem){failure.compareAndSet(null,problem);}
					}
				}
				success=failure.get()==null;
			}catch(Throwable problem){
				problem.printStackTrace(outstream);
				KillSwitch.kill("Family panel execution encountered an unrecoverable input/worker failure");
			}
		}
		private void process(ListNum<Read> batch,Scorer scorer){
			assert(batch!=null) : "A null batch is EOF and cannot contribute panel records";
			for(Read read:batch.list){
				require(read!=null && read.bases!=null && read.mate==null,"Family panels require unpaired protein sequences");
				final long before=currentCounts.queries();
				scorer.score(FamilyTailCascade.firstToken(read.id),read.bases,currentCounts);
				require(currentCounts.queries()==before+1,"Each scorer call must emit exactly one paired observation");
			}
			final long total=progress.addAndGet(batch.size()),next=nextProgress.get();
			if(total>=next && nextProgress.compareAndSet(next,total+1000000)){
				outstream.println(label+" progress: "+total+" queries scored.");
			}
		}
		final Streamer[] inputs;
		final FamilyPairedCounts[] countsT;
		FamilyPairedCounts currentCounts;
		boolean success;
	}

	private static void require(boolean valid,String why){if(!valid){throw new IllegalArgumentException(why);}}
	private static void write(Path path,ByteBuilder text){
		final ByteStreamWriter writer=new ByteStreamWriter(path.toString(),false,false,false);writer.start();
		writer.print(text);require(!writer.poisonAndWait(),"Panel summary writer failed");
	}
	private final FamilyTailCascade.LoadedClassifier classifier;
	private final Panel[] panels;
	private final int threads,inputThreads;
	private final PrintStream outstream;
	private final String label;
	final FamilyPairedCounts[] counts;
	private final ScorerFactory factory;
	private boolean errorState,started,complete;
	private final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
	private final AtomicLong progress=new AtomicLong(),nextProgress=new AtomicLong(1000000);
	private final ReadWriteLock lock=new ReentrantReadWriteLock();
	interface ScorerFactory{Scorer newScorer();}
	interface Scorer{void score(String id,byte[] raw,FamilyPairedCounts counts);}
}
