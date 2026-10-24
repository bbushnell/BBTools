package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.locks.ReadWriteLock;
import consensus.BaseGraph;
import consensus.BaseGraphHelper;
import consensus.ConsensusObject;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import idaligner.AlignmentStats;
import idaligner.QuantumAligner;
import map.ObjectIntMap;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import stream.Writer;
import stream.WriterFactory;
import structures.ByteBuilder;
import structures.IntList;
import template.Accumulator;
import template.ThreadWaiter;

/** Refines one mitochondrial consensus from every declared original member.
 * Exact duplicates share an alignment, never a vote. The seed is the largest
 * diagnostic cluster. By default broad bounds remain through stabilization;
 * roundtrim explicitly enables the secondary per-round depth-trim experiment.
 * Final trimming is followed by a new HBM in that exact coordinate frame.
 * @author Brian, Raiden */
public final class NcrnaMitoConsensus {

	public static void main(String[] args){
		final PreParser pp=new PreParser(args,NcrnaMitoConsensus.class,false);
		final Parser parser=new Parser();String seeds=null,prefix=null,name=null;
		int expected=-1,maxRounds=20;float finalTrim=0.5f,roundTrim=0;
		boolean tolerance=false;
		try{
			for(String arg:pp.args){final String[] kv=Parse.splitOnFirst(arg,'=');final String k=kv[0].toLowerCase(Locale.ROOT),v=kv.length>1?kv[1]:null;
				if(k.equals("seeds")){seeds=v;}else if(k.equals("outprefix")){prefix=v;}
				else if(k.equals("name")){name=v;}else if(k.equals("expected")){expected=Parse.parseIntKMG(v);}
				else if(k.equals("maxrounds")){maxRounds=Parse.parseIntKMG(v);}else if(k.equals("finaltrim")){finalTrim=Float.parseFloat(v);}
				else if(k.equals("roundtrim")){roundTrim=Float.parseFloat(v);}
				else if(k.equals("tolerancestable")){tolerance=Parse.parseBoolean(v);}
				else if(!parser.parse(arg,k,v)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			require(parser.in1!=null && seeds!=null && prefix!=null && name!=null && expected>0,"Require full in=, diagnostic seeds=, name=, expected= and fresh outprefix=");
			require(!parser.append && parser.maxReads<0,"All declared members must be used; read limits and append are forbidden");
			require(maxRounds>0 && Float.isFinite(finalTrim) && finalTrim>=0.49f && finalTrim<=0.51f,"Use a positive convergence cap and explicit49–51% final depth trimming");
			require(Float.isFinite(roundTrim) && (roundTrim==0 || roundTrim>=0.49f && roundTrim<=0.51f),"Explicit iteration depth fraction must be0 or49–51% before loading any input");
			require(!tolerance || roundTrim>0 && roundTrim==finalTrim,"Tolerance output must be the last per-round-trimmed sequence, with identical round/final trimming");
			RrnaResourceIO.fresh(prefix,".fa.gz",".broad.fa.gz",".hbm.gz",".rounds.tsv",".members.tsv.gz",".summary.tsv");
			Read.TO_UPPER_CASE=true;Shared.TRIM_READ_DESCRIPTION=false;ReadWrite.setZipThreads(Math.min(4,Shared.threads()));
			final ArrayList<Read> members=new ArrayList<Read>();RrnaResourceIO.stream(parser.in1,r->members.add(r));
			require(members.size()==expected,"Construction must use the full locked cohort: "+members.size()+" != "+expected);
			final ArrayList<Read> seedReads=new ArrayList<Read>();RrnaResourceIO.stream(seeds,r->seedReads.add(r));
			final Read seed=largestSeed(seedReads);final Data data=new Data(members);
			final ByteStreamWriter rounds=RrnaResourceIO.writer(prefix+".rounds.tsv");
			final Result result;
			try{rounds.println("round\treference_length\titerated_length\tcopies_voted\tunique_alignments\tstable\tcore_base_votes\treference_N_base_votes\tunanchored_copies\tconsecutive_length_change\tconsecutive_identity\ttolerance_stable");
				result=refine(data,seed.bases,name,maxRounds,finalTrim,roundTrim,Shared.threads(),rounds,tolerance);
			}finally{RrnaResourceIO.finish(rounds);}
			if(!result.exactStable){name+="_TOLERANCE_STABLE";}
			writeFasta(prefix+".fa.gz",result.consensus,name);
			writeFasta(prefix+".broad.fa.gz",result.broad,name+(roundTrim==0?"_untrimmed":"_iterated"));
			final ArrayList<BaseGraph> models=new ArrayList<BaseGraph>();models.add(result.graph);
			TrnaConsensusBuilder.writeTextModels(models,prefix+".hbm.gz");
			writeMembers(data,result.finalAlignments,prefix+".members.tsv.gz");
			final ByteStreamWriter summary=RrnaResourceIO.writer(prefix+".summary.tsv");
			try{summary.println("model\tcopies\tunique_sequences\tterminal_round\titerated_length\tfinal_length\tleft_covered_members\tright_covered_members\ttrim_rule\tround_depth_fraction\tfinal_depth_fraction\texact_stable\texact_terminal\ttolerance_stable\tconsecutive_identity\tstability");
				summary.print(new ByteBuilder().append(name).tab().append(members.size()).tab().append(data.unique.size()).tab().append(result.rounds).tab().append(result.broad.length).tab().append(result.consensus.length).tab().append(result.leftCovered).tab().append(result.rightCovered).append("\tFRACTION_OF_MAX_POSITION_DEPTH_WITH_REFERENCE_PSEUDOCOUNT").tab().append(roundTrim,3).tab().append(finalTrim,3).tab().append(result.exactStable).tab().append(result.exactTerminal).tab().append(result.toleranceStable).tab().append(result.roundIdentity,9).tab().append(result.exactStable?"EXACT_STABLE":"TOLERANCE_STABLE").nl());
			}finally{RrnaResourceIO.finish(summary);}
			pp.outstream.println("MITO_CONSENSUS_PASS copies="+members.size()+" unique="+data.unique.size()+" rounds="+result.rounds+" seed="+RrnaResourceIO.first(seed.id)+(roundTrim==0?" broad_length=":" iterated_length=")+result.broad.length+" final_length="+result.consensus.length+" votes="+result.graph.readTotal+" exact_terminal="+result.exactTerminal+" tolerance_stable="+result.toleranceStable+" stability="+(result.exactStable?"EXACT_STABLE":"TOLERANCE_STABLE")+" stable="+result.exactStable+" final_only_trim="+(roundTrim==0));
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Largest declared cluster by original copies; input order breaks equal sizes. */
	public static Read largestSeed(ArrayList<Read> seeds){
		require(seeds!=null && !seeds.isEmpty(),"A diagnostic cluster library is required, not an arbitrary input pivot");
		Read best=null;int largest=-1;
		for(Read r:seeds){final int n=Integer.parseInt(RrnaResourceIO.field(r.id,"n"));
			require(n>0 && r.bases!=null && r.length()>0,"Each diagnostic cluster must declare positive member count and sequence");
			if(n>largest){largest=n;best=r;}
		}
		return best;
	}

	/** No biological identity threshold is applied. A failed alignment is fatal. */
	public static Result refine(Data data,byte[] seed,String name,int maxRounds,float finalTrim,int threads,ByteStreamWriter rounds){
		return refine(data,seed,name,maxRounds,finalTrim,0,threads,rounds);
	}

	/** Explicit secondary experiment; zero retains the historical final-only path. */
	public static Result refine(Data data,byte[] seed,String name,int maxRounds,float finalTrim,float roundTrim,int threads,ByteStreamWriter rounds){
		return refine(data,seed,name,maxRounds,finalTrim,roundTrim,threads,rounds,false);
	}

	/** Tolerance never short-circuits the exact experiment: it is a terminal fallback. */
	public static Result refine(Data data,byte[] seed,String name,int maxRounds,float finalTrim,float roundTrim,int threads,ByteStreamWriter rounds,boolean tolerance){
		require(data!=null && seed!=null && seed.length>0 && name!=null && maxRounds>0 && threads>0,"Refinement requires a complete cohort, nonempty seed and bounded worker/iteration counts");
		require(Float.isFinite(finalTrim) && finalTrim>=0.49f && finalTrim<=0.51f,"The final consensus requires explicit49–51% endpoint trimming");
		require(Float.isFinite(roundTrim) && (roundTrim==0 || roundTrim>=0.49f && roundTrim<=0.51f),"Iteration depth trimming must be explicitly off or49–51%, never an inherited hidden cutoff");
		require(!tolerance || roundTrim>0 && roundTrim==finalTrim,"Tolerance output must retain the final iterated sequence without further depth trimming");
		final float oldTrim=ConsensusObject.trimDepthFraction,oldSub=ConsensusObject.MAF_sub,oldDel=ConsensusObject.MAF_del,oldIns=ConsensusObject.MAF_ins,oldNoRef=ConsensusObject.MAF_noref;
		final boolean oldNs=ConsensusObject.trimNs,oldMapq=ConsensusObject.useMapq,oldInvert=ConsensusObject.invertIdentity;
		try{
			ConsensusObject.trimDepthFraction=roundTrim;ConsensusObject.trimNs=true;
			ConsensusObject.MAF_sub=0.25f;ConsensusObject.MAF_del=0.5f;ConsensusObject.MAF_ins=0.5f;ConsensusObject.MAF_noref=0.4f;
			ConsensusObject.useMapq=false;ConsensusObject.invertIdentity=false;
			byte[] ref=seed.clone();final ArrayList<byte[]> history=new ArrayList<byte[]>();history.add(ref);
			for(int round=1;round<=maxRounds;round++){
				// One maximum member length on EACH side prevents a seed-length cap.
				// Native minDepth/MAF_noref may leave covered positions as N; all
				// terminal Ns are removed. Depth trimming follows explicit roundTrim.
				final BaseGraph graph=new BaseGraph(name,ref,null,0,data.maxLength);
				final Batch batch=align(data,graph.original,threads);vote(data,batch,graph);
				final byte[] broad=graph.traverse().bases;require(broad.length>0,"All-member graph produced an empty broad consensus");
				final boolean stable=Arrays.equals(ref,broad);
				final float identity=stable?1:(tolerance?consecutiveIdentity(ref,broad):Float.NaN);
				final boolean near=tolerance && withinTolerance(ref,broad,identity);
				if(rounds!=null){final ByteBuilder row=new ByteBuilder().append(round).tab().append(graph.original.length).tab().append(broad.length).tab().append(graph.readTotal).tab().append(data.unique.size()).tab().append(stable).tab().append(batch.coreVotes).tab().append(batch.nVotes).tab().append(batch.unanchoredCopies).tab().append(Math.abs(ref.length-broad.length)).tab();
					if(Float.isNaN(identity)){row.append("NA");}else{row.append(identity,9);}rounds.print(row.tab().append(near).nl());}
				boolean cycle=false;for(byte[] previous:history){cycle|=Arrays.equals(previous,broad);}
				final String terminal=stable?"EXACT":(cycle?"CYCLE":(round==maxRounds?"MAX_ROUNDS":"CONTINUE"));
				if(stable || near && (cycle || round==maxRounds)){
					// Native trimming uses a fraction of maximum position depth,
					// including reference pseudocounts, not median member endpoints.
					ConsensusObject.trimDepthFraction=finalTrim;
					final byte[] trimmed=stable?graph.traverse().bases:broad.clone();require(trimmed.length>0,"Final depth trim removed every consensus position");
					ConsensusObject.trimDepthFraction=0;
					final BaseGraph finalGraph=new BaseGraph(stable?name:name+"_TOLERANCE_STABLE",trimmed,null,0,0);
					final Batch finalBatch=align(data,trimmed,threads);vote(data,finalBatch,finalGraph);BaseGraphHelper.initForScoring(finalGraph);
					require(Arrays.equals(trimmed,finalGraph.original),"Final FASTA and HBM must have identical bases and coordinates");
					int left=0,right=0;
					for(int i=0;i<data.members.size();i++){final Read a=finalBatch.aligned[data.memberIndex.get(i)];if(a.start==0){left++;}if(a.stop==trimmed.length-1){right++;}}
					return new Result(broad,trimmed,finalGraph,finalBatch,round,left,right,stable,near,terminal,identity);
				}
				for(byte[] previous:history){require(!Arrays.equals(previous,broad),"Consensus entered a cycle at round"+round+"; do not report a nonstable resource as converged");}
				history.add(broad);ref=broad;
			}
			throw new IllegalStateException("All-member consensus did not stabilize within "+maxRounds+" rounds; inspect rounds.tsv before changing the bound");
		}finally{
			ConsensusObject.trimDepthFraction=oldTrim;ConsensusObject.trimNs=oldNs;ConsensusObject.MAF_sub=oldSub;ConsensusObject.MAF_del=oldDel;ConsensusObject.MAF_ins=oldIns;ConsensusObject.MAF_noref=oldNoRef;ConsensusObject.useMapq=oldMapq;ConsensusObject.invertIdentity=oldInvert;
		}
	}

	/** Same native shorter-query Quantum identity used in the cluster report; free ends retained. */
	public static float consecutiveIdentity(byte[] previous,byte[] current){
		require(previous!=null && current!=null && previous.length>0 && current.length>0,"Consecutive-round identity requires both complete nonempty consensus sequences");
		return NcrnaMitoSsuCompare.compare(previous,current,new AlignmentStats(true)).identity;
	}
	public static boolean withinTolerance(byte[] previous,byte[] current,float identity){
		require(previous!=null && current!=null && Float.isFinite(identity) && identity>=0 && identity<=1,"Tolerance must use a measured finite identity, never infer it from length");
		return Math.abs(previous.length-current.length)<=2 && identity>=0.999f;
	}

	private static Batch align(Data data,byte[] ref,int threads){
		assert(ref.length>0 && !data.unique.isEmpty()):"Alignment batches require the complete nonempty unique roster and current consensus";
		final Batch b=new Batch(data,ref);final int n=Math.min(threads,data.unique.size());final ArrayList<AlignThread> workers=new ArrayList<AlignThread>(n);
		for(int i=0;i<n;i++){workers.add(new AlignThread(b,i,n));}
		require(ThreadWaiter.startAndWait(workers,b),"An alignment worker failed; missing votes may never be interpreted as weak members");
		return b;
	}

	private static void vote(Data data,Batch batch,BaseGraph graph){
		assert(graph.readTotal==0 && batch.aligned.length==data.unique.size()):"Start one fresh graph from a complete alignment roster";
		for(int i=0;i<data.members.size();i++){final int u=data.memberIndex.get(i);final Read r=batch.aligned[u];require(r!=null,"Every original copy requires a successful aligned vote");graph.add(r);
			final long counts=batch.referenceCounts[u];final long core=counts>>>32;
			batch.coreVotes+=core;batch.nVotes+=counts&0xffffffffL;if(core==0){batch.unanchoredCopies++;}
		}
		require(graph.readTotal==data.members.size(),"Graph must receive exactly one vote per original copy, including seed copies and duplicates");
	}

	/** Counts diagonal query-base votes against non-N and N reference positions.
	 * High32bits=non-N core; low32bits=N. Insertions/deletions are not diagonal votes.
	 * This is telemetry only: it never excludes a member or alters its alignment. */
	public static long referenceBaseCounts(Read aligned,byte[] ref){
		require(aligned!=null && aligned.match!=null && aligned.start>=0 && aligned.stop<ref.length,"Reference telemetry requires the checked graph-frame alignment");
		int pos=aligned.start,core=0,ns=0;
		for(byte op:aligned.match){if(op=='I'){continue;}require(pos<ref.length,"Traceback escaped the reference during anchoring telemetry");
			if(op!='D'){require(op=='m' || op=='S' || op=='N',"Unexpected native traceback operation in anchoring telemetry");if(ref[pos]=='N'){ns++;}else{core++;}}pos++;
		}
		require(pos==aligned.stop+1,"Anchoring telemetry must consume the same reference span as graph voting");
		return ((long)core<<32)|(ns&0xffffffffL);
	}

	private static void writeFasta(String path,byte[] bases,String name){
		assert(bases.length>0):"A final resource cannot be empty";
		final Writer w=WriterFactory.makeWriter(FileFormat.testOutput(path,FileFormat.FASTA,null,true,false,false,true),null,1);w.start();
		final ArrayList<Read> list=new ArrayList<Read>();list.add(new Read(bases,null,name,0));w.add(list,0);
		require(!w.poisonAndWait() && w.readsWritten()==1,"Each cohort must write exactly one complete consensus");
	}

	private static void writeMembers(Data data,Batch batch,String path){
		assert(batch.aligned.length==data.unique.size()):"Member report uses final consensus alignments, not bootstrap identity";
		final ByteStreamWriter w=RrnaResourceIO.writer(path);final ByteBuilder b=new ByteBuilder();
		try{w.println("locus_id\tunique_id\tlength\tpartial\tpseudo\tquantum_identity\taligned_bases\tconsensus_start\tconsensus_stop");
			for(int i=0;i<data.members.size();i++){final Read r=data.members.get(i);final int u=data.memberIndex.get(i);final Read a=batch.aligned[u];
				w.print(b.clear().append(RrnaResourceIO.first(r.id)).tab().append(RrnaResourceIO.field(r.id,"unique")).tab().append(r.length()).tab().append(RrnaResourceIO.field(r.id,"partial")).tab().append(RrnaResourceIO.field(r.id,"pseudo")).tab().append(batch.identity[u],9).tab().append(a.length()).tab().append(a.start).tab().append(a.stop).nl());}
		}finally{RrnaResourceIO.finish(w);}
	}

	/** A retained sequence per unique ID, with a separate index for every copy. */
	public static final class Data {
		public Data(ArrayList<Read> input){
			require(input!=null && !input.isEmpty(),"An all-member consensus requires the declared nonempty cohort");members=input;
			final ObjectIntMap<String> ids=new ObjectIntMap<String>(String.class),keys=new ObjectIntMap<String>(String.class);
			for(Read r:input){require(r.bases!=null && r.length()>0 && r.mate==null,"Construction members must be nonempty unpaired sequences");
				final String id=RrnaResourceIO.first(r.id),key=RrnaResourceIO.field(r.id,"unique");require(ids.get(id)<0,"Duplicate original locus ID: "+id);ids.put(id,1);
				int u=keys.get(key);if(u<0){u=unique.size();keys.put(key,u);unique.add(r);}else{require(Arrays.equals(unique.get(u).bases,r.bases),"Duplicate-sequence key binds different bases: "+key);}
				memberIndex.add(u);maxLength=Math.max(maxLength,r.length());
			}
			require(memberIndex.size==input.size(),"Each source member needs a retained vote index");
		}
		public final ArrayList<Read> members,unique=new ArrayList<Read>();
		public final IntList memberIndex=new IntList();
		public int maxLength;
	}

	private static final class Batch implements Accumulator<AlignThread>{
		Batch(Data d,byte[] r){data=d;ref=r;aligned=new Read[d.unique.size()];identity=new float[aligned.length];referenceCounts=new long[aligned.length];}
		@Override public void accumulate(AlignThread t){ok&=t.success;}
		@Override public boolean success(){return ok;}
		@Override public ReadWriteLock rwlock(){return null;}
		final Data data;final byte[] ref;final Read[] aligned;final float[] identity;final long[] referenceCounts;boolean ok=true;
		long coreVotes,nVotes,unanchoredCopies;
	}

	private static final class AlignThread extends Thread{
		AlignThread(Batch b,int index,int stride){batch=b;first=index;step=stride;setDaemon(true);}
		@Override public void run(){
			final AlignmentStats stats=new AlignmentStats(true);
			for(int i=first;i<batch.aligned.length;i+=step){final Read member=batch.data.unique.get(i);stats.clear();stats.doTrace=true;
				final boolean swap=member.length()>batch.ref.length;
				// Quantum searches the complete reference top row. Scrabble starts a
				// narrow origin-anchored band and cannot find a seed after kilobases
				// of N padding; the large-padding regression protects this choice.
				final float identity=QuantumAligner.alignAndTraceStatic(swap?batch.ref:member.bases,swap?member.bases:batch.ref,stats);
				require(stats.matchString!=null && Float.isFinite(identity),"No usable alignment for construction member "+member.id);
				final Read aligned=swap?TrnaConsensusBuilder.invertToModelFrame(member.bases,stats):TrnaConsensusBuilder.toAlignedRead(member.bases,stats,member.id,member.numericID);
				require(aligned!=null && aligned.start>=0 && aligned.stop<batch.ref.length,"Member alignment must remain in the current graph frame: "+member.id);
				aligned.match=TrnaConsensusBuilder.normalizeOpposingIndels(aligned.match);
				require(TrnaConsensusBuilder.alignmentConsumption(aligned.match,true)==aligned.length() && TrnaConsensusBuilder.alignmentConsumption(aligned.match,false)==aligned.stop-aligned.start+1,"Traceback consumption must agree with bases and reference span before voting");
				batch.identity[i]=identity;batch.aligned[i]=aligned;batch.referenceCounts[i]=referenceBaseCounts(aligned,batch.ref);
			}
			success=true;
		}
		final Batch batch;final int first,step;boolean success;
	}

	public static final class Result{
		Result(byte[] b,byte[] c,BaseGraph g,Batch a,int n,int left,int right,boolean exact,boolean near,String terminal,float identity){broad=b;consensus=c;graph=g;finalAlignments=a;rounds=n;leftCovered=left;rightCovered=right;exactStable=exact;toleranceStable=near;exactTerminal=terminal;roundIdentity=identity;}
		public final byte[] broad,consensus;public final BaseGraph graph;public final int rounds,leftCovered,rightCovered;private final Batch finalAlignments;
		public final boolean exactStable,toleranceStable;public final String exactTerminal;public final float roundIdentity;
	}
	private static void require(boolean ok,String why){if(!ok){throw new IllegalArgumentException(why);}}
}
