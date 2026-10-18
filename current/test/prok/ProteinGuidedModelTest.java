package test.prok;

import prok.FrameStats;
import prok.GeneModel;
import prok.Orf;
import prok.ProkObject;
import prok.ProteinGuidedModel;
import prok.StatsContainer;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import dna.AminoAcid;
import dna.GeneticCode;
import stream.Read;

/** Finite independent-window oracle, frozen-copy and lifecycle tests; no HBM files or caller jobs.
 * @author Sayu
 */
public final class ProteinGuidedModelTest {

	public static void main(final String[] args) throws Exception{
		require(args.length==0, "Guided-inner fixtures accept no arguments");
		copies(); masks(); lifecycle(); partitions(); mergeRounding(); mergeLifecycle();
		System.err.println("PROTEIN_GUIDED_MODEL_TEST_PASS exact_copy=true masks=true strands=true overlaps=true "+
			"partial=true rounding=true base_totals_preserved=true immutable=true lifecycle=true overflow=true "+
			"raw_merge=true partition_bits=true merge_rounding=true merge_consumption=true merge_overflow=true");
	}

	private static void copies() throws Exception{
		final GeneModel legacy=new GeneModel(true);
		exact(legacy, legacy.copy(), false);// Null geneticCode and uncalculated -1 caches are legal snapshots.
		final GeneModel base=model(2), copy=base.copy();
		exact(base, copy, false);
		for(int c=0; c<base.allContainers.length; c++){
			final StatsContainer a=base.allContainers[c], b=copy.allContainers[c];
			require(a!=b && a.lengths!=b.lengths && a.statsArray!=b.statsArray, "Container storage must be detached");
			require(b.inner==b.statsArray[0] && b.start==b.statsArray[1] && b.stop==b.statsArray[2], "Copied aliases must index the copied FrameStats");
			for(int s=0; s<3; s++){
				final FrameStats x=a.statsArray[s], y=b.statsArray[s];
				require(x!=y && x.validSums!=y.validSums && x.counts!=y.counts && x.probs!=y.probs, "Frame storage must be detached");
				require(y.counts[0]==y.countsFalse && y.counts[1]==y.countsTrue, "Count aliases must survive copying");
				for(int f=0; f<x.frames; f++){
					require(x.countsTrue[f]!=y.countsTrue[f] && x.countsFalse[f]!=y.countsFalse[f]
						&& x.probs[f]!=y.probs[f], "Every row must be detached");
				}
			}
		}
		require(copy.baseCounts!=base.baseCounts && copy.taxIds!=base.taxIds && copy.taxIds.array!=base.taxIds.array,
			"Copied metadata must not alias cached base metadata");
		copy.statsCDS.start.probs[0][0]=99; copy.baseCounts[0]++; copy.taxIds.array[0]++;
		require(base.statsCDS.start.probs[0][0]!=99 && base.baseCounts[0]==11 && base.taxIds.get(0)==9,
			"Editing the copy must not modify its source");
		boolean rejected=false;
		try{copy.statsCDS.copyFromExact(base.statsCDS);}catch(IllegalArgumentException expected){rejected=true;}
		require(rejected, "Exact-copy destination must be fresh, not overwrite a live container");
	}

	private static void masks() throws Exception{
		final Read literal=read("literal", "ATGAAATAA");
		final List<Orf> single=Arrays.asList(orf(literal, 0, 0, 8, false, false));
		final long[][][] counts=oracle(literal, single, new boolean[]{true}, 2);
		require(sum(counts[1])==5 && sum(counts[0])==43, "Nine bases give 48 phase-windows, five positive before the real stop");
		verify(model(2), literal, single, new boolean[]{true});
		final Read seq=read("mixed", "acgtACGTNNATGAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAA");
		for(int k : new int[]{2, 3}){
			for(int phase=0; phase<3; phase++){
				for(int strand=0; strand<2; strand++){
					verify(model(k), seq, Arrays.asList(orf(seq, strand, phase, phase+14, false, false)), new boolean[]{true});
				}
			}
		}
		final Orf positive=orf(seq, 0, 12, 35, false, false);
		positive.proteinBonus=0;// Label, not bonus, defines support.
		final Orf excluded=orf(seq, 0, 12, 38, false, false);
		final Orf otherPhase=orf(seq, 0, 13, 39, false, false);
		final Orf internal=orf(seq, 0, 4, 27, true, false);
		require(internal.partialMask()==0, "Fixture must exercise internal truncation hidden by GFF edge-only partial bits");
		final Orf truncatedStop=orf(seq, 0, 24, 47, false, true);
		final Orf reverse=orf(seq, 1, 3, 29, false, false);
		final List<Orf> mixed=Arrays.asList(positive, excluded, otherPhase, internal, truncatedStop, reverse);
		final boolean[] labels={true, false, false, true, true, true};
		final GeneModel base=model(2), output=verify(base, seq, mixed, labels);
		final List<Orf> reordered=Arrays.asList(reverse, truncatedStop, internal, otherPhase, excluded, positive, positive);
		exact(output, verify(base, seq, reordered, new boolean[]{true, true, true, false, false, true, true}), false);
		// Changing only stopTruncated must protect three additional native-phase windows, not label them false.
		final Read a=read("end", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
		final Orf pos=orf(a, 0, 0, 8, false, false), end=orf(a, 0, 15, 29, false, true);
		final long[][][] extended=oracle(a, Arrays.asList(pos, end), new boolean[]{true, false}, 2);
		final Orf realEnd=orf(a, 0, 15, 29, false, false);
		final long[][][] trimmed=oracle(a, Arrays.asList(pos, realEnd), new boolean[]{true, false}, 2);
		require(sum(trimmed[0])-sum(extended[0])==3, "Truncated-stop exclusion must extend three bases past real-stop trim");
		verify(base, a, Arrays.asList(pos, end), new boolean[]{true, false});
		// Cell-first rounding: two singleton positive cells at half weight each add one, not rounded total one.
		final Read rounding=read("round", "ACGTAC");
		final GeneModel roundBase=model(2);
		final GeneModel rounded=verify(roundBase, rounding, Arrays.asList(orf(rounding, 0, 0, 5, false, false)), new boolean[]{true});
		require(rounded.statsCDS.inner.validSums[1]-roundBase.statsCDS.inner.validSums[1]==2,
			"Two singleton cells round to two; base class total must be augmented, not recomputed");
		final Read b=read("second", "ATGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCTAA");
		final List<Orf> none=Collections.emptyList();
		final ProteinGuidedModel ab=new ProteinGuidedModel(base), ba=new ProteinGuidedModel(base);
		ab.addContig(seq, mixed, labels); ab.addContig(b, none, new boolean[0]);
		ba.addContig(b, none, new boolean[0]); ba.addContig(seq, mixed, labels);
		final GeneModel both=ab.finish(); exact(both, ba.finish(), false);
		final long[][][] total=oracle(seq, mixed, labels, 2), background=oracle(b, none, new boolean[0], 2);
		for(int v=0; v<2; v++){for(int f=0; f<3; f++){for(int n=0; n<16; n++){total[v][f][n]+=background[v][f][n];}}}
		expectedInner(base, both, total, seq.length()+b.length());
		// Separate A/B/A trainers must not accumulate into the shared base.
		verify(base, b, none, new boolean[0]); exact(output, verify(base, seq, mixed, labels), false);
	}

	private static void lifecycle() throws Exception{
		final GeneModel base=model(6);
		final Read empty=read("empty", ""), shortRead=read("short", "ATGTAA"), unknown=read("unknown", "NNNNNNNNN");
		verify(base, empty, Collections.<Orf>emptyList(), new boolean[0]);
		verify(base, shortRead, Arrays.asList(orf(shortRead, 0, 0, 5, false, false)), new boolean[]{true});
		verify(base, unknown, Arrays.asList(orf(unknown, 0, 0, 8, false, false)), new boolean[]{true});
		final Read r=read("edge", "ATGAAATAA");
		for(boolean start : new boolean[]{false, true}){
			for(boolean stop : new boolean[]{false, true}){
				if(start || stop){verify(model(2), r, Arrays.asList(orf(r, 0, 0, 8, start, stop)), new boolean[]{true});}
			}
		}
		verify(model(2), r, Arrays.asList(orf(r, 0, 0, 8, false, false)), new boolean[]{false});
		final Orf rna=new Orf(r.id, 0, 8, 0, 0, r.bases, false, ProkObject.tRNA);
		verify(model(2), r, Arrays.asList(rna), new boolean[]{false});
		final ProteinGuidedModel badRna=new ProteinGuidedModel(base);
		boolean badLabel=false;
		try{badRna.addContig(r, Arrays.asList(rna), new boolean[]{true});}catch(IllegalArgumentException expected){badLabel=true;}
		require(badLabel, "A supported-stop label cannot silently label RNA as protein");
		final ProteinGuidedModel done=new ProteinGuidedModel(base);
		require(done.finish()==base, "Empty trainer must return identical base");
		boolean rejected=false;
		try{done.finish();}catch(IllegalStateException expected){rejected=true;}
		require(rejected, "finish must consume trainer exactly once");
		rejected=false; try{done.addContig(empty, Collections.<Orf>emptyList(), new boolean[0]);}catch(IllegalStateException expected){rejected=true;}
		require(rejected, "Finished trainer cannot accept another contig");
		for(int fault=0; fault<5; fault++){
			final ProteinGuidedModel trainer=new ProteinGuidedModel(base);
			final Orf o=orf(r, 0, 0, 8, false, false);
			if(fault==1){o.stop=99;} if(fault==2){o.flip();} if(fault==3){o.start=1;}
			rejected=false;
			try{trainer.addContig(fault==4 ? read("wrong", "ATGAAATAA") : r, Arrays.asList(o), fault==0 ? new boolean[0] : new boolean[]{true});}
			catch(IllegalArgumentException expected){rejected=true;}
			require(rejected, "Malformed label/coordinate/phase/contig input must fail, fault="+fault);
			rejected=false; try{trainer.finish();}catch(IllegalStateException expected){rejected=true;}
			require(rejected, "Failed trainers must not expose partial statistics");
		}
		final GeneModel overflow=model(2); overflow.statsCDS.inner.countsTrue[0][0]=Long.MAX_VALUE;
		final Read allA=read("overflow", "AAAAAAAAA");
		final ProteinGuidedModel trainer=new ProteinGuidedModel(overflow);
		trainer.addContig(allA, Arrays.asList(orf(allA, 0, 0, 8, false, false)), new boolean[]{true});
		rejected=false; try{trainer.finish();}catch(ArithmeticException expected){rejected=true;}
		require(rejected && overflow.statsCDS.inner.countsTrue[0][0]==Long.MAX_VALUE, "Overflow must fail without mutating base");
	}

	/** Vary worker assignment and reduction trees; compare raw integers before finish too. */
	private static void partitions() throws Exception{
		for(int trial=0; trial<8; trial++){
			final GeneModel base=model(2+trial%2), before=base.copy();
			final ProteinGuidedModel serial=new ProteinGuidedModel(base);
			final ProteinGuidedModel[] workers=new ProteinGuidedModel[1+trial%4];
			for(int i=0; i<workers.length; i++){workers[i]=new ProteinGuidedModel(base);}
			final Random random=new Random(314159L+trial);
			final long[][][] expected=new long[2][3][base.statsCDS.inner.kMax];
			long bases=0;
			for(int i=0; i<24; i++){
				final Read r=read("partition"+i, "acgtACGTNNATGAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAA");
				final int phase=i%3, strand=i%2;
				final List<Orf> selected=i%4==0 ? Collections.<Orf>emptyList() : Arrays.asList(
					orf(r, strand, phase, phase+26, false, false),
					orf(r, strand, phase, phase+29, false, false),
					orf(r, 1-strand, phase+3, phase+32, false, false),
					orf(r, 1-strand, phase+6, phase+35, i%2==0, true));
				final boolean[] labels=i%4==0 ? new boolean[0] : new boolean[]{true, false, false, true};
				final byte[] saved=r.bases.clone();
				serial.addContig(r, selected, labels);
				workers[random.nextInt(workers.length)].addContig(r, selected, labels);
				addRaw(expected, oracle(r, selected, labels, base.statsCDS.inner.k)); bases+=r.length();
				require(Arrays.equals(saved, r.bases), "Partition accumulation must not modify shared contig bytes");
			}
			final ArrayList<ProteinGuidedModel> tree=new ArrayList<ProteinGuidedModel>(Arrays.asList(workers));
			tree.add(new ProteinGuidedModel(base));// An idle worker still donates and becomes consumed.
			Collections.shuffle(tree, random);
			while(tree.size()>1){
				final ProteinGuidedModel donor=tree.remove(tree.size()-1), receiver=tree.get(random.nextInt(tree.size()));
				receiver.merge(donor); consumed(donor, base);
			}
			final ProteinGuidedModel merged=tree.get(0);
			rawEqual(serial, merged);
			final FrameStats raw=(FrameStats)field(merged, "delta");
			for(int v=0; v<2; v++){
				require(raw.validSums[v]==sum(expected[v]), "Raw class totals must match independent window oracle");
				for(int f=0; f<3; f++){require(Arrays.equals(raw.counts[v][f], expected[v][f]), "Raw partition counts differ from oracle");}
			}
			require((Long)field(merged, "inputBases")==bases, "Raw merge counts each input base exactly once");
			final GeneModel output=merged.finish();
			exact(serial.finish(), output, false); expectedInner(base, output, expected, bases);
			exact(before, base, false); exact(base, output, true); consumed(merged, base);
		}
	}

	private static void mergeRounding() throws Exception{
		final GeneModel base=model(2), before=base.copy();
		final Read r=read("merge-round", "ACGTAC");
		final List<Orf> selected=Arrays.asList(orf(r, 0, 0, 5, false, false));
		final boolean[] labels={true};
		final ProteinGuidedModel serial=new ProteinGuidedModel(base), receiver=new ProteinGuidedModel(base), donor=new ProteinGuidedModel(base);
		serial.addContig(r, selected, labels); serial.addContig(r, selected, labels);
		receiver.addContig(r, selected, labels); donor.addContig(r, selected, labels);
		receiver.merge(donor); rawEqual(serial, receiver);
		final GeneModel output=receiver.finish(); exact(serial.finish(), output, false);
		require(output.statsCDS.inner.validSums[1]-base.statsCDS.inner.validSums[1]==2,
			"Two worker copies of the same singleton cells must round once after merging, adding two, not four");
		consumed(donor, base);
		// A background-only worker contributes its bases to global weighting and its negative counts.
		final char[] letters=new char[240]; Arrays.fill(letters, 'A');
		final Read background=read("background", new String(letters));
		final List<Orf> none=Collections.emptyList(); final boolean[] noLabels=new boolean[0];
		final ProteinGuidedModel positive=new ProteinGuidedModel(base), negative=new ProteinGuidedModel(base), all=new ProteinGuidedModel(base);
		positive.addContig(r, selected, labels); all.addContig(r, selected, labels);
		negative.addContig(background, none, noLabels); all.addContig(background, none, noLabels);
		final FrameStats negativeBefore=((FrameStats)field(negative, "delta")).copy();
		positive.merge(negative); rawEqual(all, positive);
		// More input after a successful merge is legal, and must not alias the consumed donor's rows.
		positive.addContig(background, none, noLabels); all.addContig(background, none, noLabels);
		final FrameStats negativeAfter=(FrameStats)field(negative, "delta");
		for(int v=0; v<2; v++){for(int f=0; f<3; f++){
			require(Arrays.equals(negativeBefore.counts[v][f], negativeAfter.counts[v][f]), "Merge must not alias donor rows");
		}}
		rawEqual(all, positive);
		final long[][][] counts=oracle(r, selected, labels, 2), backgroundCounts=oracle(background, none, noLabels, 2);
		addRaw(counts, backgroundCounts); addRaw(counts, backgroundCounts);
		final GeneModel weighted=positive.finish(); exact(all.finish(), weighted, false);
		expectedInner(base, weighted, counts, r.length()+2L*background.length());
		final ProteinGuidedModel empty=new ProteinGuidedModel(base), onlyBackground=new ProteinGuidedModel(base);
		onlyBackground.addContig(background, none, noLabels); empty.merge(onlyBackground);
		empty.merge(new ProteinGuidedModel(base));
		require(empty.finish()==base, "Merging empty/background-only workers without positives must retain base identity");
		consumed(onlyBackground, base); exact(before, base, false);
	}

	private static void mergeLifecycle() throws Exception{
		final GeneModel base=model(2), before=base.copy();
		for(int fault=0; fault<4; fault++){
			final ProteinGuidedModel receiver=new ProteinGuidedModel(base);
			final ProteinGuidedModel donor=fault==0 ? null : fault==1 ? receiver
				: new ProteinGuidedModel(fault==2 ? base.copy() : base);
			if(fault==3){donor.finish();}
			boolean rejected=false;
			try{receiver.merge(donor);}catch(IllegalArgumentException | IllegalStateException expected){rejected=true;}
			require(rejected, "Null/self/different-base/consumed donor must fail, fault="+fault);
			consumed(receiver, base);
			if(fault==2){require(donor.finish()!=base, "Rejected unrelated donor must remain usable with its own base");}
		}
		// Same base identity alone is insufficient if a caller broke the immutable-geometry contract.
		final GeneModel changed=model(2);
		final ProteinGuidedModel oldGeometry=new ProteinGuidedModel(changed);
		changed.statsCDS.statsArray[0]=changed.statsCDS.inner=new FrameStats(changed.statsCDS.inner.name, 3, 3, 0);
		final ProteinGuidedModel newGeometry=new ProteinGuidedModel(changed);
		boolean rejected=false;
		try{oldGeometry.merge(newGeometry);}catch(IllegalArgumentException expected){rejected=true;}
		require(rejected, "Captured inner geometry must match even for the same base object");
		consumed(oldGeometry, changed); require(newGeometry.finish()==changed, "Geometry rejection must not consume donor");
		// Synthetic fault injection avoids enormous genomes: every checked-add layer must fail closed.
		for(int fault=0; fault<6; fault++){
			final ProteinGuidedModel receiver=new ProteinGuidedModel(base), donor=new ProteinGuidedModel(base);
			final FrameStats a=(FrameStats)field(receiver, "delta"), b=(FrameStats)field(donor, "delta");
			final long value=fault<3 ? Long.MAX_VALUE : -1;
			if(fault%3==0){
				final Field input=ProteinGuidedModel.class.getDeclaredField("inputBases"); input.setAccessible(true);
				input.setLong(receiver, value); input.setLong(donor, 1);
			}else if(fault%3==1){a.counts[1][2][15]=value; b.counts[1][2][15]=1;}
			else{a.validSums[1]=value; b.validSums[1]=1;}
			rejected=false;
			try{receiver.merge(donor);}catch(ArithmeticException | IllegalArgumentException expected){rejected=true;}
			require(rejected, "Raw base/cell/total overflow or negative value must fail, fault="+fault);
			consumed(receiver, base); consumed(donor, base); exact(before, base, false);
		}
	}

	private static void consumed(final ProteinGuidedModel trainer, final GeneModel base){
		boolean rejected=false;
		try{trainer.finish();}catch(IllegalStateException expected){rejected=true;}
		require(rejected, "Consumed trainer cannot finish");
		rejected=false;
		try{trainer.addContig(null, null, null);}catch(IllegalStateException expected){rejected=true;}
		require(rejected, "Consumed trainer cannot accumulate");
		final ProteinGuidedModel fresh=new ProteinGuidedModel(base);
		rejected=false;
		try{trainer.merge(fresh);}catch(IllegalStateException expected){rejected=true;}
		require(rejected && fresh.finish()==base, "Consumed receiver cannot merge or consume a fresh donor");
	}

	private static void rawEqual(final ProteinGuidedModel a, final ProteinGuidedModel b) throws Exception{
		require(field(a, "base")==field(b, "base") && field(a, "inputBases").equals(field(b, "inputBases")), "Raw base/length changed across partitions");
		final FrameStats x=(FrameStats)field(a, "delta"), y=(FrameStats)field(b, "delta");
		require(Arrays.equals(x.validSums, y.validSums), "Raw class totals changed across partitions");
		for(int v=0; v<2; v++){for(int f=0; f<3; f++){
			require(Arrays.equals(x.counts[v][f], y.counts[v][f]), "Raw counts changed across partitions");
		}}
	}

	private static void addRaw(final long[][][] target, final long[][][] source){
		for(int v=0; v<2; v++){for(int f=0; f<3; f++){for(int n=0; n<target[v][f].length; n++){
			target[v][f][n]=Math.addExact(target[v][f][n], source[v][f][n]);
		}}}
	}

	private static GeneModel verify(final GeneModel base, final Read r, final List<Orf> selected, final boolean[] labels) throws Exception{
		final GeneModel before=base.copy(); final byte[] saved=r.bases.clone(); final boolean[] flags=labels.clone();
		final int[] coordinates=new int[selected.size()*3];
		for(int i=0; i<selected.size(); i++){final Orf o=selected.get(i); coordinates[3*i]=o.start; coordinates[3*i+1]=o.stop; coordinates[3*i+2]=o.flipped();}
		final long[][][] counts=oracle(r, selected, labels, base.statsCDS.inner.k);
		final ProteinGuidedModel trainer=new ProteinGuidedModel(base); trainer.addContig(r, selected, labels);
		final GeneModel out=trainer.finish();
		if(sum(counts[1])==0){require(out==base, "No usable positive windows must return exact base identity, not background-only training");}
		else{require(out!=base, "Adaptation must create a detached model"); expectedInner(base, out, counts, r.length());}
		exact(before, base, false); exact(base, out, true);
		require(Arrays.equals(saved, r.bases) && Arrays.equals(flags, labels), "Training must not mutate bases or support vector");
		for(int i=0; i<selected.size(); i++){final Orf o=selected.get(i); require(coordinates[3*i]==o.start && coordinates[3*i+1]==o.stop
			&& coordinates[3*i+2]==o.flipped(), "Training must not flip shared ORFs");}
		return out;
	}

	/** Enumerate each k-window directly; no rolling state or production masks/marking helpers. */
	private static long[][][] oracle(final Read r, final List<Orf> selected, final boolean[] labels, final int k){
		final long[][][] out=new long[2][3][1<<(2*k)];
		for(int strand=0; strand<2; strand++){
			for(int end=k-1; end<r.length(); end++){
				int word=0; boolean valid=true;
				for(int p=end-k+1; p<=end; p++){
					final char c=Character.toUpperCase((char)r.bases[strand==0 ? p : r.length()-1-p]);
					int x="ACGT".indexOf(c);
					if(x<0){valid=false; break;}
					if(strand==1){x=3-x;} word=word*4+x;
				}
				if(!valid){continue;}
				for(int phase=0; phase<3; phase++){
					boolean positive=false, excluded=false;
					for(int j=0; j<selected.size(); j++){
						final Orf o=selected.get(j);
						if(o.type!=ProkObject.CDS || o.strand!=strand){continue;}
						final int first=strand==0 ? o.start : r.length()-1-o.stop;
						final int last=strand==0 ? o.stop : r.length()-1-o.start;
						if(end-k+1<first || end>last-(o.stopTruncated ? 0 : 3) || (end-first)%3!=phase){continue;}
						if(labels[j] && !o.startTruncated && !o.stopTruncated && o.length()%3==0){positive=true;}
						else{excluded=true;}
					}
					if(positive || !excluded){out[positive ? 1 : 0][phase][word]++;}
				}
			}
		}
		return out;
	}

	private static void expectedInner(final GeneModel base, final GeneModel actual, final long[][][] raw, final long bases){
		final FrameStats b=base.statsCDS.inner, a=actual.statsCDS.inner;
		final double weight=0.50*Math.max(100L, base.basesProcessed)/Math.max(100L, bases);
		final long[] sums=b.validSums.clone();
		for(int v=0; v<2; v++){for(int f=0; f<3; f++){for(int n=0; n<b.kMax; n++){
			final long addition=Math.round(raw[v][f][n]*weight); sums[v]+=addition;
			require(a.counts[v][f][n]==b.counts[v][f][n]+addition, "Window oracle differs at class/phase/kmer "+v+"/"+f+"/"+n);
		}}}
		require(Arrays.equals(sums, a.validSums), "New cell-first sums must add to unchanged base totals");
		final float average=(float)((sums[1]+1.0)/(sums[0]+sums[1]+1.0)), inv=1.0f/average;
		for(int f=0; f<3; f++){for(int n=0; n<b.kMax; n++){
			final float prob=(float)(a.countsTrue[f][n]/(a.countsTrue[f][n]+a.countsFalse[f][n]+1.0))*inv;
			bits(prob, a.probs[f][n], "Updated inner probabilities must be calculated after merging");
		}}
	}

	/** Deliberately noncanonical cached values detect accidental recalculation of frozen blocks. */
	private static GeneModel model(final int k) throws Exception{
		final GeneModel m=new GeneModel(false, GeneticCode.forTable(4));
		for(StatsContainer sc : m.allContainers){
			sc.setInner(sc.type==ProkObject.CDS ? k : 1, sc.type==ProkObject.CDS ? 3 : 1);
			sc.setStart(1, 2, 1); sc.setStop(1, 2, 0);
			for(FrameStats fs : sc.statsArray){
				for(int v=0; v<2; v++){for(int f=0; f<fs.frames; f++){for(int n=0; n<fs.kMax; n++){
					fs.counts[v][f][n]=v+f+n%3; fs.validSums[v]+=v+f+n%3;
				}}}
				fs.validSums[0]+=7; fs.validSums[1]+=11;// Existing totals are NOT silently repaired.
				fs.calculate(); fs.probs[0][0]=-0f; fs.probs[0][1]=Float.intBitsToFloat(0x7fc12345);
			}
			sc.lengths[5]=17; sc.lengthSum=123; sc.lengthCount=7; sc.lengthAvg=91; sc.invLengthAvg=-0f;
		}
		m.numFiles=3; m.readsProcessed=4; m.basesProcessed=100; m.genesProcessed=5; m.filesProcessed=6;
		m.baseCounts=new long[]{11, 12, 13, 14, 15}; m.taxIds.add(9); m.taxIds.add(2);
		final Field max=GeneModel.class.getDeclaredField("maxReads"); max.setAccessible(true); max.setLong(m, 1234);
		return m;
	}

	private static void exact(final GeneModel a, final GeneModel b, final boolean skipInner) throws Exception{
		require(a.numFiles==b.numFiles && a.readsProcessed==b.readsProcessed && a.basesProcessed==b.basesProcessed
			&& a.genesProcessed==b.genesProcessed && a.filesProcessed==b.filesProcessed, "Model metadata changed");
		require(Arrays.equals(a.baseCounts, b.baseCounts) && a.taxIds.size==b.taxIds.size, "Base/taxonomy metadata changed");
		for(int i=0; i<a.taxIds.size; i++){require(a.taxIds.get(i)==b.taxIds.get(i), "Taxonomy ordering changed");}
		for(String name : new String[]{"maxReads", "geneticCode"}){
			final Object x=field(a, name), y=field(b, name);
			require(x==y || (x!=null && x.equals(y)), "Private model metadata changed: "+name);
		}
		for(int c=0; c<a.allContainers.length; c++){
			final StatsContainer x=a.allContainers[c], y=b.allContainers[c];
			require(x.type==y.type && x.name.equals(y.name) && x.lengthSum==y.lengthSum && x.lengthCount==y.lengthCount
				&& x.lengthAvg==y.lengthAvg && Arrays.equals(x.lengths, y.lengths), "Length/type statistics changed");
			bits(x.invLengthAvg, y.invLengthAvg, "Frozen inverse length changed");
			for(int s=0; s<3; s++){
				final FrameStats p=x.statsArray[s], q=y.statsArray[s];
				require(p.compatibleWith(q), "Frame geometry changed");
				if(skipInner && c==0 && s==0){continue;}
				require(Arrays.equals(p.validSums, q.validSums), "Frozen totals changed");
				bits((Float)field(p, "average"), (Float)field(q, "average"), "Cached average changed");
				bits((Float)field(p, "invAvg"), (Float)field(q, "invAvg"), "Cached inverse average changed");
				for(int f=0; f<p.frames; f++){
					require(Arrays.equals(p.countsTrue[f], q.countsTrue[f]) && Arrays.equals(p.countsFalse[f], q.countsFalse[f]), "Frozen count row changed");
					for(int n=0; n<p.kMax; n++){bits(p.probs[f][n], q.probs[f][n], "Frozen probability changed");}
				}
			}
		}
	}
	private static Object field(final Object o, final String name) throws Exception{
		final Field f=o.getClass().getDeclaredField(name); f.setAccessible(true); return f.get(o);
	}
	/** Keep literal lowercase/empty fixture bytes, independent of global Read normalization flags. */
	private static Read read(final String id, final String sequence){return new Read(sequence.getBytes(StandardCharsets.US_ASCII), null, id, 0, false);}
	private static Orf orf(final Read r, final int strand, final int first, final int last, final boolean start, final boolean stop){
		return new Orf(r.id, first, last, strand, first%3, strand==0 ? r.bases : AminoAcid.reverseComplementBases(r.bases), true, ProkObject.CDS, start, stop);
	}
	private static long sum(final long[][] a){long sum=0; for(long[] row : a){for(long n : row){sum+=n;}} return sum;}
	private static void bits(final float a, final float b, final String why){require(Float.floatToRawIntBits(a)==Float.floatToRawIntBits(b), why);}
	private static void require(final boolean ok, final String why){if(!ok){throw new AssertionError(why);}}
}
