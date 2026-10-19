package test.prok;

import prok.*;

import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;
import dna.AminoAcid;
import stream.SamLine;

/** Real full-model comparisons and actual GeneCaller DP, including an intentional
 * neither-SSU fixture that must block verification rather than trigger a rollback.
 * @author Raiden */
public final class SsuOverlapResolverTest {
	public static void main(String[] args) throws Exception{
		final int oldLimit=SamLine.INTRON_LIMIT;
		final ArrayList<NcrnaFamily> oldFamilies=new ArrayList<NcrnaFamily>(GeneCaller.ncrnaFamilies);
		try{
			SamLine.INTRON_LIMIT=Integer.MAX_VALUE;
			GeneCaller.ncrnaFamilies.clear();
			GeneCaller.ncrnaFamilies.add(new NcrnaFamily("euk18S",new byte[][]{random(300,101)},null,
				new String[]{"euk18S_mitochondria"},null,17,100,100));
			for(int strand=0; strand<2; strand++){partialAndPath(strand);}
			tiesAndEligibility();betterMito();tandem();cacheReset();neitherBlocks();alternateSsuBlocks();
		}finally{SamLine.INTRON_LIMIT=oldLimit;GeneCaller.ncrnaFamilies.clear();GeneCaller.ncrnaFamilies.addAll(oldFamilies);}
		System.out.println("SSU_VETO_TEST_PASS full_models shared_span both_strands partial tie no16S primary_mito tandem cache final_DP neither_blocks alternate_SSU_blocks");
	}

	static void partialAndPath(int strand) throws Exception{
		final byte[] forward=random(800,31),model=Arrays.copyOfRange(forward,100,400),mito=extend(model,300);
		final byte[] bases=strand==0?forward:AminoAcid.reverseComplementBases(forward);
		final int start=strand==0?100:400,stop=start+299;
		final Orf s=ssu(bases,start,stop,strand,model),m=mito(bases,start,stop,strand,mito,true);
		final ArrayList<Orf>[] lists=lists(s,m);
		check(path(lists,bases).contains(m),"Higher original mito score must win the actual baseline DP fixture");
		final SsuOverlapResolver r=new SsuOverlapResolver(false,true);
		r.apply(lists,bases);
		check(lists[strand].size()==1 && lists[strand].get(0)==s,"A clipped perfect local match must not beat the better full16S model");
		check(s.start==start && s.stop==stop && s.orfScore==500f,"Veto must preserve16S endpoints and global score scale");
		final ArrayList<Orf> selected=path(lists,bases);
		check(selected.contains(s) && !selected.contains(m),"Preferred16S must actually survive global DP");
		r.finish(selected);check(r.preferredRetained==1 && r.neitherRetained==0 && r.alignmentCount==2,"Both actual full-model alignments and final retention must be counted");
		check(Arrays.equals(bases,strand==0?forward:AminoAcid.reverseComplementBases(forward)),"Comparison cannot change genomic bases");
	}

	static void tiesAndEligibility() throws Exception{
		final byte[] bases=random(800,91),model=Arrays.copyOfRange(bases,100,400);
		for(int mode=0; mode<4; mode++){
			final Orf s=ssu(bases,100,399,0,model),m=mito(bases,100,399,mode==2?1:0,model.clone(),mode!=3);
			final ArrayList<Orf>[] lists=mode==1?lists(m):lists(s,m);
			final SsuOverlapResolver r=new SsuOverlapResolver(false,true);r.apply(lists,bases);
			check(lists[m.strand].contains(m),"Tie/no16S/opposite-strand/primary-mito cases must preserve the mitochondrial candidate");
			check(r.vetoCount==0 && r.alignmentCount==(mode==0?2:0),"Ineligible pairs must not align; exact identity ties must not veto");
			r.finish(path(lists,bases));
		}
	}

	static void betterMito() throws Exception{
		final byte[] bases=random(800,19),model=Arrays.copyOfRange(bases,100,400),weaker=model.clone();
		for(int i=5; i<weaker.length; i+=4){weaker[i]=weaker[i]=='A'?(byte)'C':(byte)'A';}
		final Orf s=ssu(bases,100,399,0,weaker),m=mito(bases,100,399,0,model,true);
		final ArrayList<Orf>[] lists=lists(s,m);final SsuOverlapResolver r=new SsuOverlapResolver(false,true);r.apply(lists,bases);
		check(lists[0].contains(m) && r.vetoCount==0,"The better full mitochondrial model must remain eligible despite a competing16S");
		check(path(lists,bases).contains(m),"Mitochondrial control must still reach final output");r.finish(path(lists,bases));
	}

	static void tandem(){
		final byte[] bases=random(800,43),model=Arrays.copyOfRange(bases,0,300),bad=random(500,71);
		final Orf a=ssu(bases,0,299,0,model),b=mito(bases,150,449,0,bad,true),c=mito(bases,300,599,0,bad,true);
		check(SsuOverlapResolver.overlap(a,b) && SsuOverlapResolver.overlap(b,c) && !SsuOverlapResolver.overlap(a,c),"Fixture must exercise a direct-half-overlap chain");
		final ArrayList<Orf>[] lists=lists(a,b,c);final SsuOverlapResolver r=new SsuOverlapResolver(false,false);r.apply(lists,bases);
		check(!lists[0].contains(b) && lists[0].contains(c),"Veto must not propagate transitively through an adjacent mitochondrial candidate");
		r.finish(new ArrayList<Orf>(Arrays.asList(a,c)));
	}

	static void cacheReset(){
		final byte[] a=random(300,199),b=random(300,411),bases=random(800,213);
		System.arraycopy(a,0,bases,100,300);
		final SsuOverlapResolver r=new SsuOverlapResolver(false,true);
		final Orf s=ssu(bases,100,399,0,a),m=mito(bases,100,399,0,b,true),copy=mito(bases,100,399,0,b,true);
		final ArrayList<Orf>[] first=lists(m,s,copy);r.apply(first,bases);
		check(r.alignmentCount==2 && r.vetoCount==2 && first[0].size()==1,"Repeated model/span fits must reuse cached alignments without duplicating16S candidates");r.finish(first[0]);
		System.arraycopy(b,0,bases,100,300);
		final Orf s2=ssu(bases,100,399,0,a),m2=mito(bases,100,399,0,b,true);final ArrayList<Orf>[] next=lists(s2,m2);r.apply(next,bases);
		check(r.alignmentCount==4 && r.vetoCount==2 && next[0].contains(m2),"Cached coordinates from another contig must not reverse the new full-model winner");r.finish(next[0]);
	}

	static void neitherBlocks() throws Exception{
		final byte[] bases=random(800,55),model=Arrays.copyOfRange(bases,100,400);
		final Orf s=ssu(bases,100,399,0,model),m=mito(bases,100,399,0,extend(model,300),true);
		final ArrayList<Orf>[] lists=lists(s,m);final SsuOverlapResolver r=new SsuOverlapResolver(false,true);r.apply(lists,bases);
		final Orf cds=new Orf("fixture",100,399,0,0,bases,false,ProkObject.CDS);cds.orfScore=2000f;lists[0].add(cds);
		final ArrayList<Orf> selected=path(lists,bases);check(selected.contains(cds) && !selected.contains(s),"Fixture must make the preferred16S lose the actual global DP");
		boolean blocked=false;try{r.finish(selected);}catch(IllegalStateException e){blocked=e.getMessage().contains("Neither SSU");}
		check(blocked && r.neitherRetained==1,"A real neither-SSU outcome must block validation without silently restoring mito");
	}

	static void alternateSsuBlocks() throws Exception{
		final byte[] bases=random(800,76),model=Arrays.copyOfRange(bases,100,400),weaker=model.clone();
		for(int i=7; i<weaker.length; i+=6){weaker[i]=weaker[i]=='A'?(byte)'C':(byte)'A';}
		final Orf preferred=ssu(bases,100,399,0,model),alternative=ssu(bases,100,399,0,weaker),m=mito(bases,100,399,0,extend(model,300),true);
		alternative.orfScore=700f;
		final ArrayList<Orf>[] lists=lists(preferred,alternative,m);final SsuOverlapResolver r=new SsuOverlapResolver(false,true);r.apply(lists,bases);
		final ArrayList<Orf> selected=path(lists,bases);check(selected.contains(alternative) && !selected.contains(preferred),"Fixture must retain an alternative SSU while losing the preferred full-model fit");
		boolean blocked=false;try{r.finish(selected);}catch(IllegalStateException e){blocked=e.getMessage().contains("Preferred 16S");}
		check(blocked && r.neitherRetained==0 && r.preferredLost==1,"Another SSU surviving cannot masquerade as retention of the preferred16S");
	}

	@SuppressWarnings("unchecked")
	static ArrayList<Orf> path(ArrayList<Orf>[] rnas, byte[] bases) throws Exception{
		final ArrayList<Orf>[] frames=new ArrayList[6];for(int i=0; i<6; i++){frames[i]=new ArrayList<Orf>();}
		for(ArrayList<Orf> list:rnas){for(Orf o:list){o.prevPlus=o.prevMinus=null;o.pathScorePlus=o.pathScoreMinus=-999999;frames[3*o.strand+o.frame].add(o);}}
		for(ArrayList<Orf> list:frames){Collections.sort(list);}
		final GeneCaller caller=new GeneCaller(60,1,1,0f,0f,0f,0f,0f,new GeneModel(false));
		final Method method=GeneCaller.class.getDeclaredMethod("findPath",ArrayList[].class,byte[].class);method.setAccessible(true);
		return (ArrayList<Orf>)method.invoke(caller,frames,bases);
	}
	@SuppressWarnings("unchecked")
	static ArrayList<Orf>[] lists(Orf... orfs){final ArrayList<Orf>[] result=new ArrayList[]{new ArrayList<Orf>(),new ArrayList<Orf>()};for(Orf o:orfs){result[o.strand].add(o);}return result;}
	static Orf ssu(byte[] bases,int start,int stop,int strand,byte[] model){final Orf o=new Orf("fixture",start,stop,strand,0,bases,false,ProkObject.r16S);o.orfScore=500f;o.ssuFullModel=model;o.ssuFullModelName="fixture16S";return o;}
	static Orf mito(byte[] bases,int start,int stop,int strand,byte[] model,boolean rescue){final Orf o=new Orf("fixture",start,stop,strand,0,bases,false,ProkObject.r18S);o.orfScore=900f;o.ncrnaFamily="euk18S";o.ssuFullModel=model;o.ssuFullModelName="euk18S_mitochondria";o.ssuClippedMitoRescue=rescue;return o;}
	static byte[] random(int length,int seed){final Random r=new Random(seed);final byte[] out=new byte[length],alphabet={'A','C','G','T'};for(int i=0; i<length; i++){out[i]=alphabet[r.nextInt(4)];}return out;}
	static byte[] extend(byte[] sequence,int tail){final byte[] out=Arrays.copyOf(sequence,sequence.length+tail),extra=random(tail,989);System.arraycopy(extra,0,out,sequence.length,tail);return out;}
	static void check(boolean ok,String why){if(!ok){throw new AssertionError(why);}}
}
