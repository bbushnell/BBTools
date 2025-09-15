package var;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import dna.Data;
import fileIO.TextStreamWriter;
import shared.Parse;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;

public class StackVariations2 {
	
	public static void main(String[] args){
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			//outstream=pp.outstream;
		}
		
		Timer t=new Timer();
		
		String inPattern=(args[0].equalsIgnoreCase("null") ? null : args[0]);
		String outPattern=args[1];
		
		assert(!inPattern.equalsIgnoreCase(outPattern));
		
		int minChrom=-1;
		int maxChrom=-1;
		
		boolean filter=false;
		Data.GENOME_BUILD=-1;
		
		for(int i=2; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equalsIgnoreCase("filter")){
				filter=true;
			}else if(a.startsWith("filter")){
				if(b.equals("1") || b.startsWith("t")){filter=true;}
				else if(b.equals("0") || b.startsWith("f")){filter=false;}
				else{throw new RuntimeException("Unknown parameter "+args[i]);}
			}else if(a.equalsIgnoreCase("strict")){
				if(b==null){STRICT=true;}
				else if(b.equals("1") || b.startsWith("t")){STRICT=true;}
				else if(b.equals("0") || b.startsWith("f")){STRICT=false;}
				else{throw new RuntimeException("Unknown parameter "+args[i]);}
			}else if(a.equals("genome") || a.equals("build")){
				Data.setGenome(Integer.parseInt(b));
				if(minChrom==-1){minChrom=1;}
				if(maxChrom==-1){maxChrom=Data.numChroms;}
			}else if(a.equals("minchrom")){
				minChrom=Integer.parseInt(b);
			}else if(a.equals("maxchrom")){
				maxChrom=Integer.parseInt(b);
			}else if(a.equals("threads") || a.equals("t")){
				THREADS=Integer.parseInt(b);
			}else if(a.equals("minreads")){
				MIN_READS_TO_KEEP=Integer.parseInt(b);
			}else if(a.equals("blocksize")){
				GenerateVarlets2.BLOCKSIZE=(Integer.parseInt(b));
			}else if(a.equals("deletefiles") || a.startsWith("deletetemp") || a.startsWith("deleteinput") || a.equals("delete")){
				DELETE_INPUT=(Parse.parseBoolean(b));
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		assert(minChrom>=0 && maxChrom>=minChrom) : "Please set minchrom and maxchrom.";
		if(Data.GENOME_BUILD<0){throw new RuntimeException("Please set genome number.");}
		THREADS=Tools.max(1, THREADS);
		
//		for(byte i=minChrom; i<=maxChrom; i++){
//			String fname1=inPattern.replace("#", i+"");
//			String fname2=outPattern.replace("#", i+"");
//			assert(new File(fname1).exists());
//			assert(!new File(fname2).exists());
//			processFile(fname1, fname2, filter);
//		}
		
		runThreaded(inPattern, outPattern, minChrom, maxChrom, filter);
		
		t.stop();
		System.out.println("Input Vars:        \t"+(totalIn_global-totalInNR_global));
		System.out.println("Input No-ref:      \t"+totalInNR_global);
		System.out.println("Input Delta Length:\t"+deltaLenIn_global);
		System.out.println();
		System.out.println("Kept Vars:         \t"+(totalKept_global-totalKeptNR_global));
		System.out.println("Kept No-ref:       \t"+totalKeptNR_global);
		System.out.println("Kept Snp:          \t"+snpKept_global);
		System.out.println("Kept Del:          \t"+delKept_global+"\t\tLength: \t"+delLenKept_global);
		System.out.println("Kept Ins:          \t"+insKept_global+"\t\tLength: \t"+insLenKept_global);
		System.out.println("Kept Sub:          \t"+subKept_global+"\t\tLength: \t"+subLenKept_global);
		System.out.println("Kept Delta Length: \t"+deltaLenKept_global);
		System.out.println("Kept Avg Score:    \t"+(scoreKept_global/(Tools.max(1, totalKept_global))));
		System.out.println();
		System.out.println("Dropped Vars:      \t"+(totalDropped_global-totalDroppedNR_global));
		System.out.println("Dropped No-ref:    \t"+totalDroppedNR_global);
		System.out.println("Dropped Avg Score: \t"+(scoreDropped_global/Tools.max(1, totalDropped_global)));
		System.out.println();
		System.out.println("Time: \t"+t);
	}
	
	public static final void runThreaded(String inPattern, String outPattern, int minChrom, int maxChrom, boolean filter){
		ArrayList<SVThread> svts=new ArrayList<SVThread>();
		for(int i=minChrom; i<=maxChrom; i++){
			assert(inPattern==null || !inPattern.equalsIgnoreCase(outPattern));
			String fname1=inPattern;
			String fname2=outPattern.replace("#", i+"");
			addThread(1);
			SVThread svt=new SVThread(fname1, fname2, i, filter);
			svts.add(svt);
			new Thread(svt).start();
		}
		while(addThread(0)>0){}
		for(SVThread svt : svts){

			snpKept_global+=svt.snpKept;
			delKept_global+=svt.delKept;
			insKept_global+=svt.insKept;
			subKept_global+=svt.subKept;
			delLenKept_global+=svt.delLenKept;
			insLenKept_global+=svt.insLenKept;
			subLenKept_global+=svt.subLenKept;
			deltaLenKept_global+=svt.deltaLenKept;

			deltaLenIn_global+=svt.deltaLenIn;
			totalIn_global+=svt.totalIn;
			totalInNR_global+=svt.totalInNR;
			totalKept_global+=svt.totalKept;
			totalDropped_global+=svt.totalDropped;
			totalKeptNR_global+=svt.totalKeptNR;
			totalDroppedNR_global+=svt.totalDroppedNR;
			scoreKept_global+=svt.scoreKept;
			scoreDropped_global+=svt.scoreDropped;
		}
	}
	
	
	public static boolean passesFilterSNP(Varlet v){
		

			//Best so far:

			if(STRICT){
				
				if(v.endDist<3){return false;}
				if(v.tailDist<10){return false;}

				//NOTE!  Last thing I did was make this more strict by adding 1 to all the num reads/unique reads required.
				if(v.minStrandReads()>=2){
					
					if(v.errors>2){return false;}
					if(v.expectedErrors>1.5f){return false;}
//					if(v.expectedErrors-v.errors>3f){return false;}
					if(v.maxReadQuality()<18){return false;}
					if(v.avgReadQuality()<13){return false;}
					if(v.maxVarQuality()<26){return false;}
					if(v.avgVarQuality()<18){return false;}
					if(v.numReads<4){return false;}
					if(v.numSemiUniqueReads<4){return false;}
					if(v.numUniqueReads<2){return false;}
					if(v.paired<3){return false;}
					
				}else if(v.minStrandReads()>=1){
					
					if(v.errors>2){return false;}
					if(v.expectedErrors>1.2f){return false;}
//					if(v.expectedErrors-v.errors>3f){return false;}
					if(v.maxReadQuality()<19){return false;}
					if(v.avgReadQuality()<14){return false;}
					if(v.maxVarQuality()<28){return false;}
					if(v.avgVarQuality()<19){return false;}
					if(v.numReads<3){return false;}
					if(v.numSemiUniqueReads<3){return false;}
					if(v.numUniqueReads<2){return false;}
					if(v.paired<3){return false;}
					
				}else{
					if(v.endDist<8){return false;}
					if(v.tailDist<14){return false;}
					
					if(v.errors>0){return false;}
					if(v.expectedErrors>0.5f){return false;}
//					if(v.expectedErrors-v.errors>2f){return false;}
					if(v.maxReadQuality()<21){return false;}
					if(v.avgReadQuality()<17){return false;}
					if(v.maxVarQuality()<30){return false;}
					if(v.avgVarQuality()<21){return false;}
					if(v.numReads<6){return false;}
					if(v.numSemiUniqueReads<5){return false;}
					if(v.numUniqueReads<3){return false;}
					if(v.paired<5){return false;}
					if(v.score()<8100){return false;}
				}
				
//				else{
//					if(v.endDist<8){return false;}
//					if(v.tailDist<14){return false;}
//
//					if(v.errors>0){return false;}
//					if(v.expectedErrors>0.5f){return false;}
////					if(v.expectedErrors-v.errors>2f){return false;}
//					if(v.maxReadQuality()<21){return false;}
//					if(v.avgReadQuality()<17){return false;}
//					if(v.maxVarQuality()<30){return false;}
//					if(v.avgVarQuality()<21){return false;}
//					if(v.numReads<5){return false;}
//					if(v.numSemiUniqueReads<4){return false;}
//					if(v.numUniqueReads<2){return false;}
//					if(v.paired<4){return false;}
//					if(v.score()<8100){return false;}
//				}

			}else{
				
				assert(false) : "disabled";

			}
			
		
		
		return true;
	}
	
	public static boolean passesFilterOther(Varlet v){
		
		
			
				if(v.endDist<3){return false;}
				if(v.tailDist<10){return false;}

				//NOTE!  Last thing I did was make this more strict by adding 1 to all the num reads/unique reads required.
				if(v.minStrandReads()>=2){
					
					if(v.errors>2){return false;}
					if(v.expectedErrors>1.5f){return false;}
//					if(v.expectedErrors-v.errors>3f){return false;}
					if(v.maxReadQuality()<16){return false;}
					if(v.avgReadQuality()<12){return false;}
					if(v.maxVarQuality()<26){return false;}
					if(v.avgVarQuality()<16){return false;}
					if(v.numReads<4){return false;}
					if(v.numSemiUniqueReads<4){return false;}
					if(v.numUniqueReads<2){return false;}
					if(v.paired<3){return false;}
					
				}else if(v.minStrandReads()>=1){
					
					if(v.errors>2){return false;}
					if(v.expectedErrors>1.2f){return false;}
//					if(v.expectedErrors-v.errors>3f){return false;}
					if(v.maxReadQuality()<17){return false;}
					if(v.avgReadQuality()<13){return false;}
					if(v.maxVarQuality()<28){return false;}
					if(v.avgVarQuality()<17){return false;}
					if(v.numReads<4){return false;}
					if(v.numSemiUniqueReads<4){return false;}
					if(v.numUniqueReads<2){return false;}
					if(v.paired<3){return false;}
					
				}else{
					if(v.endDist<8){return false;}
					if(v.tailDist<14){return false;}
					
					if(v.errors>0){return false;}
					if(v.expectedErrors>0.5f){return false;}
//					if(v.expectedErrors-v.errors>2f){return false;}
					if(v.maxReadQuality()<20){return false;}
					if(v.avgReadQuality()<16){return false;}
					if(v.maxVarQuality()<30){return false;}
					if(v.avgVarQuality()<20){return false;}
					if(v.numReads<6){return false;}
					if(v.numSemiUniqueReads<5){return false;}
					if(v.numUniqueReads<3){return false;}
					if(v.paired<5){return false;}
					if(v.score()<6500){return false;}
				}
				
			
			
		
		
		return true;
	}
	
	
	public static ArrayList<Varlet> mergeAll(ArrayList<Varlet> vars){
		if(vars==null || vars.size()==0){return null;}
		ArrayList<Varlet> out=new ArrayList<Varlet>(8+vars.size()/16);
		Shared.sort(vars);
		
		ArrayList<Varlet> temp=new ArrayList<Varlet>(64);
		for(int i=0; i<vars.size(); i++){
//			while(vars.get(i).beginLoc<3746582){i++;}
			Varlet v=vars.get(i);
//			System.err.println("Grabbed "+v.beginLoc+" ~ "+v.call);
			if(temp.isEmpty()){
//				System.err.println("Adding "+v.beginLoc+" ~ "+v.call);
				temp.add(v);
			}else{
				if(v.equals(temp.get(0))){
					temp.add(v);
//					System.err.println("Adding "+v.beginLoc+" ~ "+v.call);
				}else{
//					System.err.println("Merging "+temp.size()+" x "+v.beginLoc+" ~ "+v.call);
					Varlet result=mergeEqualVarlets(temp);
					if(result.numReads>MIN_READS_TO_KEEP){
						out.add(result);
					}else if(result.numReads==MIN_READS_TO_KEEP){
						if(result.maxVarQuality()>=MIN_QUALITY_AT_MIN_READS &&
								result.errors<=MAX_ERRORS_AT_MIN_READS &&
								result.expectedErrors<=MAX_EXPECTED_ERRORS_AT_MIN_READS &&
								(result.paired>0 || !REQUIRE_PAIRED_AT_MIN_READS)){
							out.add(result);
						}
					}
					temp.clear();
					temp.add(v);
				}
			}
			
			
		}
		
		if(!temp.isEmpty()){
			if(temp.size()>=MIN_READS_TO_KEEP){
				Varlet result=mergeEqualVarlets(temp);
				out.add(result);
			}
			temp.clear();
		}
		
		{//For testing
			Shared.sort(out); //Should already be sorted...
			for(int i=1; i<out.size(); i++){
				assert(!out.get(i).equals(out.get(i-1)));
			}
		}
		
		
		if(verbose){System.err.println("out.size="+out.size());}
		
		return out;
	}
	
	
	public static Varlet mergeEqualVarlets(ArrayList<Varlet> vars){
		
//		System.err.println("Merging "+vars.size()+" vars.");
		
		if(vars.size()==1){return vars.get(0);}

		HashMap<Integer, ArrayList<Varlet>> plus=new HashMap<Integer, ArrayList<Varlet>>(Tools.min(8, vars.size()));
		HashMap<Integer, ArrayList<Varlet>> minus=new HashMap<Integer, ArrayList<Varlet>>(Tools.min(8, vars.size()));
		
		int numReads=0;
		int numSemiUniqueReads=0;
		int numUniqueReads=0;
		int pairedReads=0;
		int plusReads1=0;
		int minusReads1=0;
		int plusReads2=0;
		int minusReads2=0;

		int totalQuality=0;
		int totalVarQuality=0;
		
		int maxReadQuality=0;
		int maxVarQuality=0;

		int maxMapScore=0;
		int bestLen=0;
		int minReadStart=Integer.MAX_VALUE;
		int maxReadStop=-999999;
		
		int maxHeadDist=-1;
		int maxTailDist=-1;
		int maxEndDist=-1;
		
		Varlet bestVar=null;
		
		int minErrors=999;
		float minExpectedErrors=999f;
		
		for(Varlet v : vars){
			
			numReads+=v.numReads;
			numSemiUniqueReads+=v.numSemiUniqueReads;
			plusReads1+=v.numPlusReads1;
			minusReads1+=v.numMinusReads1;
			plusReads2+=v.numPlusReads2;
			minusReads2+=v.numMinusReads2;
			
			if(v.errors<minErrors || (v.errors<=minErrors && v.maxReadQuality()>maxReadQuality)){
				bestVar=v;
			}
			
			totalQuality+=v.avgReadQuality()*v.numReads;
			maxReadQuality=Tools.max(maxReadQuality, v.maxReadQuality());
			
			totalVarQuality+=v.avgVarQuality()*v.numReads;
			maxVarQuality=Tools.max(maxVarQuality, v.maxVarQuality());
			
			if(bestLen==0 || (v.mapScore>=maxMapScore && v.readLen>=bestLen)){
				bestLen=v.readLen;
			}

			maxHeadDist=Tools.max(maxHeadDist, v.headDist);
			maxTailDist=Tools.max(maxTailDist, v.tailDist);
			maxEndDist=Tools.max(maxEndDist, v.endDist);
			
			minErrors=Tools.min(minErrors, v.errors);
			minExpectedErrors=Tools.min(minExpectedErrors, v.expectedErrors);
			maxMapScore=Tools.max(maxMapScore, v.mapScore);
			minReadStart=Tools.min(minReadStart, v.readStart);
			maxReadStop=Tools.max(maxReadStop, v.readStop);
			assert(minReadStart<maxReadStop) : "\n"+minReadStart+"\n"+maxReadStop+"\n"+v.toText();
			
			pairedReads+=v.paired;
			
			if(v.strand==Shared.PLUS){
				ArrayList<Varlet> value=plus.get(v.readStart);
				if(value==null){
					numUniqueReads++;
					value=new ArrayList<Varlet>(2);
					plus.put(v.readStart, value);
				}
				value.add(v);
			}else{
				ArrayList<Varlet> value=minus.get(v.readStop);
				if(value==null){
					numUniqueReads++;
					value=new ArrayList<Varlet>(2);
					minus.put(v.readStop, value);
				}
				value.add(v);
			}
		}
		
//		byte plusReads=(byte) ((plus.isEmpty() ? 0 : 1)+(minus.isEmpty() ? 0 : 1));
		
		float avgVarQuality=totalVarQuality/(float)numReads;
		float avgReadQuality=totalQuality/(float)numReads;

		int netQuality=(int)Math.ceil((avgVarQuality+maxVarQuality)/2);
		int netReadQuality=(int)Math.ceil((avgReadQuality+maxReadQuality)/2);
		
		Varlet v=new Varlet(bestVar.chromosome, ((plusReads1+plusReads2>0) && (minusReads1+minusReads2>0) ? Shared.PLUS : bestVar.strand),
				bestVar.beginLoc, bestVar.endLoc, bestVar.matchStart, bestVar.matchStop, bestVar.varType, bestVar.ref, bestVar.call,
				netQuality, netReadQuality, maxMapScore, minErrors, minExpectedErrors, pairedReads, bestVar.readID, bestLen,
				minReadStart, maxReadStop, numReads, maxHeadDist, maxTailDist, maxEndDist, bestVar.pairNum());
		
		
		v.setMaxReadQuality(maxReadQuality);
		v.setMaxVarQuality(maxVarQuality);
		v.setAvgReadQuality((int)Math.ceil(avgReadQuality));
		v.setAvgVarQuality((int)Math.ceil(avgVarQuality));
		
		v.numSemiUniqueReads=numSemiUniqueReads;
		v.numUniqueReads=numUniqueReads;
		v.numPlusReads1=plusReads1;
		v.numMinusReads1=minusReads1;
		v.numPlusReads2=plusReads2;
		v.numMinusReads2=minusReads2;
		assert(plusReads1+minusReads1+plusReads2+minusReads2==numSemiUniqueReads);
		
		assert(v.numReads>=v.numSemiUniqueReads);
		assert(v.numSemiUniqueReads>=v.numUniqueReads);
		
		//This assertion is only correct if stacking is done from raw, uncombined varlets.
		assert(v.numSemiUniqueReads==vars.size()) : "\n"+vars.size()+", "+v.numReads+", "+v.numSemiUniqueReads+", "+v.numUniqueReads
			+"\n"+v.toText();
		
		assert(v.numUniqueReads<=v.numReads && v.numUniqueReads>0);
		assert(v.numUniqueReads==plus.size()+minus.size()) : "numUniqueReads="+numUniqueReads+
		", v.numUniqueReads="+v.numUniqueReads+", v.numReads="+v.numReads
		+", plus.size()="+plus.size()+", minus.size()="+minus.size()+"\n"+vars+"\n";
		
		return v;
	}
	
	
	private static class SVThread implements Runnable {
		
		public SVThread(String fname1_, String fname2_, final int chrom_, boolean filter_){
			fname1=fname1_;
			fname2=fname2_;
			filter=filter_;
			chrom=chrom_;
		}
		
		@Override
		public void run() {
//			addThread(1);
			assert(activeThreads>0);
			processFile(fname1, fname2);
			addThread(-1);
		}
		
		private final void processFile(final String inName, final String outName){
			
			final long[] keys=GenerateVarlets2.keys(chrom);
			final TextStreamWriter tsw=(inName==null ? null : new TextStreamWriter(outName, true, false, false));
			if(tsw!=null){
				tsw.start();
				tsw.println(Varlet.textHeader());
			}
			
			for(final long key : keys){
				String blockname=GenerateVarlets2.fname(key, inName);
				
				ArrayList<Varlet> initial=Varlet.fromTextFile(blockname);

				for(Varlet v : initial){
					if(v.varType==Variation.NOREF){totalInNR++;}
					totalIn++;

					int dif=v.lengthDif();
					deltaLenIn+=dif;
				}

				if(verbose){System.err.println("Initial:  \t"+initial.size());}
				
				int merged=mergeAll2(initial, tsw);
				
				initial=null;
				if(verbose){System.err.println("Merged:   \t"+merged);}
				
			}
			
			if(tsw!=null){
				tsw.poison();
				if(DELETE_INPUT){
					for(int i=0; i<10 && tsw.isAlive(); i++){
						try {
							tsw.join(10000);
						} catch (InterruptedException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
					if(tsw.isAlive()){
						System.err.println(tsw.getClass().getName()+" for "+outName+" refused to die.");
						assert(false);
					}
				}
			}
			
			if(DELETE_INPUT){
				for(final long key : keys){
					String blockname=GenerateVarlets2.fname(key, inName);
//					System.out.println("Deleting "+blockname);
					new File(blockname).delete();
				}
			}
		}
		
		
		
		
		
		private final int mergeAll2(ArrayList<Varlet> vars, TextStreamWriter tsw){
			if(vars==null || vars.size()==0){return 0;}
			
			Shared.sort(vars);
			int out=0;
			
			ArrayList<Varlet> temp=new ArrayList<Varlet>(64);
			for(int i=0; i<vars.size(); i++){
//				while(vars.get(i).beginLoc<3746582){i++;}
//				Varlet v=vars.get(i);
				final Varlet v=vars.set(i, null);
//				System.err.println("Grabbed "+v.beginLoc+" ~ "+v.call);
				if(temp.isEmpty()){
//					System.err.println("Adding "+v.beginLoc+" ~ "+v.call);
					temp.add(v);
				}else{
					if(v.equals(temp.get(0))){
						temp.add(v);
//						System.err.println("Adding "+v.beginLoc+" ~ "+v.call);
					}else{
//						System.err.println("Merging "+temp.size()+" x "+v.beginLoc+" ~ "+v.call);
						Varlet result=mergeEqualVarlets(temp);
						
						processMergedVar(result, tsw);
						out++;
						
						temp.clear();
						temp.add(v);
					}
				}
			}
			
			if(!temp.isEmpty()){
				if(temp.size()>=MIN_READS_TO_KEEP){
					Varlet result=mergeEqualVarlets(temp);
					out++;
					processMergedVar(result, tsw);
				}
				temp.clear();
			}
			
			return out;
		}
		
		
		private final boolean processMergedVar(Varlet v, TextStreamWriter tsw){
			
			if(v==null){return false;}
			if(v.numReads<MIN_READS_TO_KEEP){return false;}
			if(v.numReads==MIN_READS_TO_KEEP){
				if(v.maxVarQuality()<MIN_QUALITY_AT_MIN_READS ||
						v.errors<=MAX_ERRORS_AT_MIN_READS ||
						v.expectedErrors<=MAX_EXPECTED_ERRORS_AT_MIN_READS ||
						(v.paired<1 && REQUIRE_PAIRED_AT_MIN_READS)){
					return false;
				}
			}
			
			boolean keep;
			
			if(filter){
				keep=filterLight(v);
			}else{
				keep=true;
				totalKept++;
				scoreKept+=v.score();
			}
			
			if(keep){
				StringBuilder sb=v.toText();
				sb.append('\n');
				tsw.print(sb);
			}
			return keep;
		}
		
		
		private final boolean filterLight(Varlet v){
			int dropped=0;

			int dif=v.lengthDif();
//			deltaLenIn+=dif;

			boolean passes=true;
			if(v.varType==Variation.NOCALL){
				passes=false;
			}else if(v.numSemiUniqueReads<2){
				passes=false;
			}else if(v.endDist<6 || v.tailDist<10){
				passes=false;
			}else if(v.maxVarQuality()<24){
				passes=false;
			}else if(v.expectedErrors>2){
				passes=false;
			}

			if(passes && STRICT){
				passes=passesFilterLight(v);
			}

			if(passes){
				if(v.varType==Variation.NOREF){totalKeptNR++;}
				else if(v.varType==Variation.SNP){snpKept++;}
				else if(v.varType==Variation.DEL){
					delKept++;
					//						delLenKept-=v.lengthRef();
					delLenKept+=dif;
				}
				else if(v.varType==Variation.INS){
					insKept++;
					//						insLenKept+=v.lengthVar();
					insLenKept+=dif;
				}
				else if(v.varType==Variation.DELINS){
					subKept++;
					//						subLenKept+=(v.lengthRef()-v.lengthVar());
					subLenKept+=dif;
				}
				totalKept++;
				scoreKept+=v.score();
				deltaLenKept+=dif;
			}else{
				if(v.varType==Variation.NOREF){totalDroppedNR++;}
				dropped++;
				scoreDropped+=v.score();
			}

			totalDropped+=dropped;
			return passes;
		}
		
		private static boolean passesFilterLight(Varlet v){
			if(v.endDist<4){return false;}
			if(v.tailDist<10){return false;}
			
			//NOTE!  Last thing I did was make this more strict by adding 1 to all the num reads/unique reads required.
			if(v.minStrandReads()>=2){

				if(v.errors>2){return false;}
				if(v.expectedErrors>1.4f){return false;}
				//			if(v.expectedErrors-v.errors>3f){return false;}
				if(v.maxReadQuality()<17){return false;}
				if(v.avgReadQuality()<13){return false;}
				if(v.maxVarQuality()<26){return false;}
				if(v.avgVarQuality()<17){return false;}
				if(v.numReads<3){return false;}
				if(v.numSemiUniqueReads<3){return false;}
				if(v.numUniqueReads<2){return false;}
//				if(v.paired<3){return false;}
				if(v.score()<8200){return false;}

			}else if(v.minStrandReads()>=1){
				if(v.endDist<7){return false;}
				if(v.tailDist<12){return false;}

				if(v.errors>2){return false;}
				if(v.expectedErrors>1.1f){return false;}
				//			if(v.expectedErrors-v.errors>3f){return false;}
				if(v.maxReadQuality()<18){return false;}
				if(v.avgReadQuality()<14){return false;}
				if(v.maxVarQuality()<28){return false;}
				if(v.avgVarQuality()<18){return false;}
				if(v.numReads<4){return false;}
				if(v.numSemiUniqueReads<3){return false;}
				if(v.numUniqueReads<2){return false;}
//				if(v.paired<3){return false;}
				if(v.score()<8020){return false;}
			}else{
				if(v.endDist<8){return false;}
				if(v.tailDist<14){return false;}

				if(v.errors>0){return false;}
				if(v.expectedErrors>0.5f){return false;}
				//			if(v.expectedErrors-v.errors>2f){return false;}
				if(v.maxReadQuality()<21){return false;}
				if(v.avgReadQuality()<17){return false;}
				if(v.maxVarQuality()<30){return false;}
				if(v.avgVarQuality()<21){return false;}
				if(v.numReads<6){return false;}
				if(v.numSemiUniqueReads<5){return false;}
				if(v.numUniqueReads<3){return false;}
//				if(v.paired<5){return false;}
				if(v.score()<7670){return false;}
			}
			return true;
		}
		
		private long deltaLenKept=0;
		private long snpKept=0;
		private long delKept=0;
		private long insKept=0;
		private long subKept=0;
		private long delLenKept=0;
		private long insLenKept=0;
		private long subLenKept=0;
		
		private long deltaLenIn=0;
		private long totalIn=0;
		private long totalInNR=0;
		
		private long totalKept=0;
		private long totalKeptNR=0;
		private long totalDropped=0;
		private long totalDroppedNR=0;
		private long scoreKept=0;
		private long scoreDropped=0;
		
		private final String fname1;
		private final String fname2;
		private final boolean filter;
		private final int chrom;
	}
	
	private static int addThread(int x){
		synchronized(THREADLOCK){
			while(x>0 && activeThreads>=THREADS){
				try {
					THREADLOCK.wait(200);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			activeThreads+=x;
			return activeThreads;
		}
	}


	public static long deltaLenKept_global=0;
	public static long deltaLenIn_global=0;
	
	public static long snpKept_global=0;
	public static long delKept_global=0;
	public static long insKept_global=0;
	public static long subKept_global=0;
	public static long delLenKept_global=0;
	public static long insLenKept_global=0;
	public static long subLenKept_global=0;
	
	public static long totalIn_global=0;
	public static long totalInNR_global=0;
	public static long totalKept_global=0;
	public static long totalDropped_global=0;
	public static long totalKeptNR_global=0;
	public static long totalDroppedNR_global=0;
	public static long scoreKept_global=0;
	public static long scoreDropped_global=0;
	
	private static int activeThreads=0;
	
	private static final String THREADLOCK=new String("THREADLOCK");
	private static int THREADS=7;
	private static boolean DELETE_INPUT=false;
	public static int MIN_READS_TO_KEEP=1;
	public static final int MIN_QUALITY_AT_MIN_READS=14;
	public static final int MAX_ERRORS_AT_MIN_READS=2;
	public static final int MAX_EXPECTED_ERRORS_AT_MIN_READS=4;
	public static final boolean REQUIRE_PAIRED_AT_MIN_READS=false;
	public static boolean STRICT=false;
	public static boolean VSTRICT=false;
	public static boolean USTRICT=false;
	
	public static final boolean verbose=false;
}
