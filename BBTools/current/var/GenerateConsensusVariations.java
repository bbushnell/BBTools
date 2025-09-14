package var;

import dna.ChromosomeArray;
import dna.Data;
import fileIO.ReadWrite;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.PreParser;
import shared.Timer;
import shared.Tools;
import structures.CoverageArray;

/**
 * @author Brian Bushnell
 * @date Jul 23, 2012
 *
 */
public class GenerateConsensusVariations {
	
	public static void main(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		Timer t=new Timer();

		String inVarsPattern=args[0];
		String inCovPattern=args[1];
		String outPattern=args[2];

		assert(!inVarsPattern.equalsIgnoreCase(outPattern));
		assert(!inCovPattern.equalsIgnoreCase(outPattern));
		
		int minChrom=-1;
		int maxChrom=-1;
		int minCoverage=1;
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		
		for(int i=3; i<args.length; i++){
			final String arg=args[i].toLowerCase();
			String[] split=arg.split("=");
			String a=split[0];
			String b=split.length>1 ? split[1] : null;
			
			if(a.startsWith("mincov")){
				minCoverage=Integer.parseInt(b);
				assert(minCoverage>0);
			}else if(a.startsWith("consensus")){
				consensusRatio=Float.parseFloat(b);
//				assert(consensusRatio>=0.5f && consensusRatio<=1f);
				assert(consensusRatio>=0f && consensusRatio<=1f);
				consensusRatioNR=1-(1-consensusRatio)*.5f; //Lower multiplier is more accurate
//				assert(false) : consensusRatioNR;
			}else if(a.equals("genome") || a.equals("build")){
				Data.setGenome(Integer.parseInt(b));
				if(minChrom==-1){minChrom=1;}
				if(maxChrom==-1){maxChrom=Data.numChroms;}
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
			}else if(a.equals("minchrom")){
				minChrom=Integer.parseInt(b);
			}else if(a.equals("maxchrom")){
				maxChrom=Integer.parseInt(b);
			}else if(a.equals("threads") || a.equals("t")){
				THREADS=Integer.parseInt(b);
			}else if(a.startsWith("noref") || a.startsWith("undef")){
				NOREF_CAP=Integer.parseInt(b);
			}else{
				System.err.println("Unknown argument "+arg);
			}
		}
		
		for(int chrom=minChrom; chrom<=maxChrom; chrom++){
			process(inVarsPattern.replaceFirst("#", ""+chrom), inCovPattern.replaceFirst("#", ""+chrom), outPattern.replaceFirst("#", ""+chrom), chrom, minCoverage);
		}
		
		t.stop();
		
		System.out.println();
		System.out.println("Vars in:          \t"+(VARS_IN-NOREFS_IN));
		System.out.println("Length Delta in:  \t"+VARLEN_IN);
		System.out.println("No-refs in:       \t"+NOREFS_IN);
		System.out.println();
		System.out.println("Vars out:         \t"+(VARS_OUT-NOREFS_OUT));
		System.out.println("Length Delta out: \t"+VARLEN_OUT);
		System.out.println("No-refs out:      \t"+NOREFS_OUT);
		System.out.println();
		System.out.println("Time: \t"+t);
		
	}
	
	/** Now removes overlapping vars by retaining better quality one. */
	public static void process(final String invars, final String incov, final String outfile, final int chrom, final int mincov){
		TextFile tf=new TextFile(invars, true);
		CoverageArray ca=ReadWrite.read(CoverageArray.class, incov, true);
		TextStreamWriter tsw=new TextStreamWriter(outfile, true, false, true);
		tsw.start();
		
		ChromosomeArray cha=Data.getChromosome(chrom);
		
		Varlet prev=null;
		
		tsw.println(Varlet.header());
		
		for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
			if(s.charAt(0)!='#'){
				Varlet v=Varlet.fromText(s);
				VARS_IN++;
				int dif=v.lengthDif();
				VARLEN_IN+=dif;
				if(v.varType==Variation.NOREF){NOREFS_IN++;}
				
				boolean passes=passesFilter(v, ca, cha, mincov);
				boolean overlap=(prev==null ? false : v.beginLoc<=prev.endLoc);
//				if(passes){System.out.println(v.varTypeMap[v.varType]+"  " +
//						((v.ref==null || v.ref.length()<1 ? "." : v.ref)+" "+(v.call==null || v.call.length()<1 ? "." : v.call))+
//						"  \tchr"+v.chromosome+" "+v.beginLoc+"        \tdepth "+v.numReads+" / "+ca.get(v.beginLoc)+"");}
				
				if(!overlap){
					if(prev!=null){
						StringBuilder sb=prev.toText().append('\n');
						tsw.print(sb);
						VARS_OUT++;
						VARLEN_OUT+=prev.lengthDif();
						if(prev.varType==Variation.NOREF){NOREFS_OUT++;}
					}
					prev=null;
				}else{
					if(passes && v.score()>prev.score()){
						prev=null;
					}else{
						v=null;
					}
				}
				
				if(passes && v!=null){
					prev=v;
				}
				
//				if(passesFilter(v, ca, cha, mincov)){
//					StringBuilder sb=v.toText().append('\n');
//					tsw.print(sb);
//					VARS_OUT++;
//					VARLEN_OUT+=dif;
//					if(v.varType==Variation.NOREF){NOREFS_OUT++;}
//				}else{
//
//				}
			}
		}
		
		if(prev!=null){
			StringBuilder sb=prev.toText().append('\n');
			tsw.print(sb);
			VARS_OUT++;
			VARLEN_OUT+=prev.lengthDif();
			if(prev.varType==Variation.NOREF){NOREFS_OUT++;}
		}
		
		tf.close();
		tsw.poison();
		Data.unload(chrom, true);
		
	}
	
	
	/**
	 * @param v
	 * @param ca
	 * @return
	 */
	private static boolean passesFilter(Varlet v, CoverageArray ca, ChromosomeArray cha, int minCoverageToPass) {
		
		int dif=v.lengthDif();
		
		int midLoc=(v.beginLoc+v.endLoc)/2;
		int midCov=ca.get(midLoc);
		int maxCov=midCov, minCov=midCov;
		
		int bound1, bound2;
		float ratio;
		
		if(verbose){System.err.println("\nConsidering varlet "+v);}
		
		if(v.varType==Variation.NOREF){
			bound1=v.beginLoc;
			bound2=v.endLoc;
			minCoverageToPass=minCoverageToPass*2+5;
			ratio=consensusRatioNR;
		}else{
			bound1=v.beginLoc;
			bound2=v.endLoc;
			ratio=consensusRatio;
//			if(dif<0){minCoverageToPass++;} //Helps reduce deletion bias
		}
		
		for(int i=bound1; i<=bound2; i++){
			int cov=ca.get(i);
			minCov=Tools.min(minCov, cov);
			maxCov=Tools.max(maxCov, cov);
			if(verbose){System.err.println("minCov = "+minCov+", maxCov = "+maxCov);}
		}
//		if(dif<)
		
		if(minCov<minCoverageToPass){
			if(verbose){System.err.println("Low coverage, "+minCov+"<"+minCoverageToPass+"\n"+v);}
			return false;
		}
		int minReads=(int)Math.ceil(ratio*minCov);
		if(v.numReads<minReads){
			if(verbose){System.err.println("Low reads, mincov="+minCov+", "+v.numReads+"<"+minReads+"\n"+v);}
			return false;
		}
		if(v.minStrandReads()<1 && v.numSemiUniqueReads<2*minCoverageToPass){
			if(verbose){System.err.println("Low strands, mincov="+minCov+", "+v.minStrandReads()+"<"+1+"\n"+v);}
			return false;
		}
		
		//Check noref
		if(v.varType==Variation.NOREF){
			if(NOREF_CAP>=0){
				int a=Tools.max(v.beginLoc-NOREF_CAP, cha.minIndex);
				int b=Tools.min(v.endLoc+NOREF_CAP, cha.maxIndex);
				if(cha.isFullyUndefined(a, b)){
					if(verbose){System.err.println("Noref cap, mincov="+minCov+"\n"+v);}
					return false;
				}
			}
		}
		if(verbose){System.err.println("Retaining variation.");}
		return true;
	}


	/** TODO */
	public static int THREADS=1;
	public static int NOREF_CAP=-1;
	public static float consensusRatio=1f;
	public static float consensusRatioNR=1f;
	public static long VARS_IN=0;
	public static long VARLEN_IN=0;
	public static long NOREFS_IN=0;
	public static long VARS_OUT=0;
	public static long VARLEN_OUT=0;
	public static long NOREFS_OUT=0;
	public static boolean verbose=false;
	
}
