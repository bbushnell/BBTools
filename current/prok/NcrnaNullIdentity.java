package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.ReadWriteLock;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import idaligner.QuantumAligner;
import idaligner.ScrabbleAligner;
import map.ObjectSet;
import parse.PreParser;
import shared.Shared;
import stream.Read;
import structures.ByteBuilder;
import template.Accumulator;
import template.ThreadWaiter;

/** Seed-bypassed null identity distribution, never a production caller FPR.
 * Uses NcrnaScavenger's initial generic verifier direction and length dispatch:
 * Quantum(model, sequence) above quantum; Scrabble(sequence, model) otherwise.
 * Both orientations and every model are measured without seed, index, span or HBM gates.
 * @author Raiden
 */
public final class NcrnaNullIdentity implements Accumulator<NcrnaNullIdentity.Worker> {
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, NcrnaNullIdentity.class, false);String input=null, ref=null, out=null, label=null;
		int expected=-1, expectedModels=-1, threads=1, quantum=-1;
		for(String arg:pp.args){
			final int eq=arg.indexOf('=');require(eq>0, "Use flag=value");final String k=arg.substring(0, eq).toLowerCase(Locale.ROOT), v=arg.substring(eq+1);
			if(k.equals("in")){input=v;}else if(k.equals("ref")){ref=v;}else if(k.equals("outprefix")){out=v;}else if(k.equals("label")){label=v;}
			else if(k.equals("expected")){expected=Integer.parseInt(v);}else if(k.equals("expectedmodels")){expectedModels=Integer.parseInt(v);}
			else if(k.equals("t")){threads=Integer.parseInt(v);}else if(k.equals("quantum")){quantum=Integer.parseInt(v);}else{throw new IllegalArgumentException("Unknown option: "+k);}
		}
		require(input!=null && ref!=null && out!=null && label!=null && expected>0 && expectedModels>0 && threads>0 && quantum>0,
			"Bind exact sample, complete library, denominators, explicit caller dispatch threshold and interpretation label");
		require(label.equals("NULL_PANEL") || label.equals("NULL_NONPANEL") || label.equals("POSITIVE_CONTROL"), "Keep null and positive interpretation explicit");
		Shared.setThreads(threads);Shared.TRIM_READ_DESCRIPTION=false;Read.TO_UPPER_CASE=true;
		try{new NcrnaNullIdentity(input, ref, expected, expectedModels, quantum).process(out, threads, label);}finally{Shared.closeStream(pp.outstream);}
	}
	private NcrnaNullIdentity(String input, String ref, int expected, int expectedModels, int q){
		quantum=q;library=TrnaConsensusBuilder.loadLibrary(ref);require(library!=null && library.length==expectedModels, "Every declared consensus must participate; no shortlist cap");
		final ObjectSet<String> ids=new ObjectSet<String>(String.class);
		RrnaResourceIO.stream(input, r->{
			final String id=RrnaResourceIO.first(r.id);require(r.bases!=null && r.length()>0 && ids.add(id), "Null distribution requires unique nonempty records");
			for(byte base:r.bases){if(RrnaResourceIO.base(base)<0){throw new IllegalArgumentException("Ambiguous core outside the declared sampling policy: "+id);}}
			members.add(new Member(id, r.bases));
		});
		require(members.size()==expected, "Complete sample denominator changed: "+members.size()+" expected="+expected);
	}
	private void process(String prefix, int threads, String label){
		RrnaResourceIO.fresh(prefix, ".loci.tsv", ".quantiles.tsv", ".thresholds.tsv", ".cm_review.ids.txt", ".summary.tsv");
		final ArrayList<Worker> workers=new ArrayList<Worker>();for(int i=0; i<Math.min(threads, members.size()); i++){workers.add(new Worker(this));}
		if(!ThreadWaiter.startAndWait(workers, this)){throw new IllegalStateException("Null identity worker failed; no results published", failure);}
		require(calls==2L*library.length*members.size(), "Both orientations of every member require exactly one alignment to every model");
		write(prefix, label);
	}
	@Override public void accumulate(Worker w){assert(!w.isAlive()):"ThreadWaiter must join workers before accumulating scores";calls+=w.calls;failed|=!w.success;if(failure==null){failure=w.failure;}}
	@Override public boolean success(){return !failed;}
	@Override public ReadWriteLock rwlock(){return null;}
	private void write(String prefix, String label){
		final ByteStreamWriter loci=RrnaResourceIO.writer(prefix+".loci.tsv"), quantiles=RrnaResourceIO.writer(prefix+".quantiles.tsv"), thresholds=RrnaResourceIO.writer(prefix+".thresholds.tsv"), review=RrnaResourceIO.writer(prefix+".cm_review.ids.txt"), summary=RrnaResourceIO.writer(prefix+".summary.tsv");
		final ByteBuilder b=new ByteBuilder();final float[][] values=new float[3][members.size()];final long[][] above=new long[3][CUTS.length];long bases=0;
		try{
			loci.println("id\tlabel\tlength\tbest_forward\tbest_reverse\tbest_either\twinning_strand\tmodel_index\taligned_start0\taligned_stop0\tinterpretation");
			for(int i=0; i<members.size(); i++){
				final Member m=members.get(i);bases+=m.bases.length;final int strand=m.best[1]>m.best[0]?1:0;final float best=m.best[strand];
				values[0][i]=m.best[0];values[1][i]=m.best[1];values[2][i]=best;
				for(int view=0; view<3; view++){for(int c=0; c<CUTS.length; c++){if(values[view][i]>=CUTS[c]){above[view][c]++;}}}
				loci.print(b.clear().append(m.id).tab().append(label).tab().append(m.bases.length).tab().appendSlow(m.best[0]).tab().appendSlow(m.best[1]).tab().appendSlow(best)
					.tab().append(strand==0?"+":"-").tab().append(m.model[strand]).tab().append(m.start[strand]).tab().append(m.stop[strand]).tab().append(INTERPRETATION).nl());
				if(best>=.60f){review.println(m.id);}
			}
			quantiles.println("label\torientation\tmembers\tmin\tp10\tmedian\tp90\tp99\tmax\tquantile_rule\tinterpretation");
			thresholds.println("label\torientation\tidentity_threshold\tmembers\treaching_threshold\tfraction\tinterpretation");
			for(int view=0; view<3; view++){
				final float[] sorted=values[view];Arrays.sort(sorted);
				b.clear().append(label).tab().append(VIEWS[view]).tab().append(sorted.length);
				for(double q:new double[]{0, .1, .5, .9, .99, 1}){b.tab().appendSlow(sorted[(int)Math.floor((sorted.length-1)*q)]);}
				quantiles.print(b.append("\tFLOOR_P_TIMES_N_MINUS_ONE\t").append(INTERPRETATION).nl());
				for(int c=0; c<CUTS.length; c++){
					thresholds.print(b.clear().append(label).tab().append(VIEWS[view]).tab().appendSlow(CUTS[c]).tab().append(sorted.length).tab().append(above[view][c])
						.tab().append(above[view][c]/(double)sorted.length, 8).tab().append(INTERPRETATION).nl());
				}
			}
			summary.println("label\tmembers\tbases\tmodels\torientations\talignments\tquantum_length_threshold\tinterpretation");
			summary.print(b.clear().append(label).tab().append(members.size()).tab().append(bases).tab().append(library.length).append("\t2\t").append(calls).tab().append(quantum).tab().append(INTERPRETATION).nl());
		}finally{RrnaResourceIO.finish(loci, quantiles, thresholds, review, summary);}
		System.err.println("NCRNA_NULL_IDENTITY_PASS members="+members.size()+" models="+library.length+" alignments="+calls+"; "+INTERPRETATION);
	}
	static final class Worker extends Thread {
		Worker(NcrnaNullIdentity x){owner=x;}
		@Override public void run(){
			try{
				for(int first; (first=owner.next.getAndAdd(8))<owner.members.size(); ){
					for(int i=first, end=Math.min(first+8, owner.members.size()); i<end; i++){score(owner.members.get(i));}
				}success=true;
			}catch(Throwable t){failure=t;owner.next.set(owner.members.size());}
		}
		private void score(Member m){
			final byte[] reverse=m.bases.clone();AminoAcid.reverseComplementBasesInPlace(reverse);
			for(int strand=0; strand<2; strand++){
				final byte[] seq=strand==0?m.bases:reverse;float best=-1;
				for(int model=0; model<owner.library.length; model++){
					final float identity;
					if(seq.length>owner.quantum){identity=QuantumAligner.alignStatic(owner.library[model], seq, pos);}
					else{identity=ScrabbleAligner.alignStatic(seq, owner.library[model], null);pos[0]=0;pos[1]=seq.length-1;}
					require(Float.isFinite(identity) && identity>=0 && identity<=1, "Native initial-verifier identity must be finite in [0,1]");calls++;
					if(identity>best){best=identity;m.model[strand]=model;m.start[strand]=pos[0];m.stop[strand]=pos[1];}
				}m.best[strand]=best;
			}
		}
		final NcrnaNullIdentity owner;final int[] pos=new int[4];long calls;boolean success;Throwable failure;
	}
	private static final class Member {
		Member(String i, byte[] b){id=i;bases=b;}final String id;final byte[] bases;final float[] best=new float[2];final int[] model=new int[2], start=new int[2], stop=new int[2];
	}
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
	private final byte[][] library;private final int quantum;private final ArrayList<Member> members=new ArrayList<Member>();private final AtomicInteger next=new AtomicInteger();
	private long calls;private boolean failed;private Throwable failure;
	private static final float[] CUTS={.60f, .61f, .62f, .63f, .64f, .65f};
	private static final String[] VIEWS={"FORWARD", "REVERSE", "BEST_OF_BOTH"};
	public static final String INTERPRETATION="NULL_IDENTITY_DISTRIBUTION_SEED_BYPASSED_NOT_CALLER_FPR";
}
