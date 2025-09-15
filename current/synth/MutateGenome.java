package synth;

import java.util.ArrayList;
import java.util.Random;

import bin.ConservationModel;
import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import shared.Parse;
import shared.Parser;
import shared.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ConcurrentReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.FASTQ;
import stream.Read;
import structures.ByteBuilder;
import structures.ListNum;
import var2.Var;

/**
 * @author Brian Bushnell
 * @date June 1, 2016
 *
 */
public class MutateGenome {

	public static void main(String[] args){
		if(parsePloidy(args)!=1){
			MutateGenome2.main(args);
			return;
		}
		Timer t=new Timer();
		MutateGenome x=new MutateGenome(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public MutateGenome(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		Shared.setBufferLen(1);
		FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		
		Parser parser=new Parser();
		parser.overwrite=true;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("period")){
				period=Integer.parseInt(b);
			}else if(a.equals("subrate") || a.equals("snprate")){
				subRate=Float.parseFloat(b);
				if(subRate>1){subRate/=100;}
			}else if(a.equals("indelrate")){
				float indelRate=Float.parseFloat(b);
				if(indelRate>1){indelRate/=100;}
				insRate=delRate=indelRate/2;
			}else if(a.equals("insrate")){
				insRate=Float.parseFloat(b);
				if(insRate>1){insRate/=100;}
			}else if(a.equals("delrate")){
				delRate=Float.parseFloat(b);
				if(delRate>1){delRate/=100;}
			}else if(a.equals("maxindel")){
				maxIndel=Parse.parseIntKMG(b);
			}else if(a.equals("indelspacing")){
				indelSpacing=Parse.parseIntKMG(b);
			}else if(a.equals("seed")){
				seed=Long.parseLong(b);
			}else if(a.equals("wavecov") || a.equals("sinewave")){
				sinewaves=(Parse.parseBoolean(b) ? 3 : 0);
			}else if(a.equals("waves") || a.equals("sinewaves")){
				sinewaves=Integer.parseInt(b);
			}else if(a.equals("pad")){
				padLeft=padRight=Integer.parseInt(b);
			}else if(a.equals("padleft")){
				padLeft=Integer.parseInt(b);
			}else if(a.equals("padright")){
				padRight=Integer.parseInt(b);
			}else if(a.equals("nohomopolymers") || a.equals("nohomop") || 
					a.equals("banhomopolymers") || a.equals("banhomop")){
				banHomopolymers=Parse.parseBoolean(b);
			}else if(a.equals("prefix")){
				prefix=b;
			}else if(a.equals("vcf") || a.equals("outvcf") || a.equals("vcfout")
					|| a.equals("vars") || a.equals("varsout") || a.equals("outvars")){
				outVcf=b;
			}else if(a.equals("id") || a.equals("identity")){
				float x=Float.parseFloat(b);
				if(x>1){x=x/100;}
				x=1-x;
				insRate=delRate=x*0.005f;
				subRate=x-insRate-delRate;
			}else if(a.equals("fraction") || a.equals("completeness")){
				float x=Float.parseFloat(b);
				if(x>1){x=x/100;}
				genomeFraction=x;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		float indelRate=insRate+delRate;
		errorRate=subRate+indelRate;
		errorRate2=subRate+(indelRate*Math.max(1f, 0.16666f*(maxIndel+1)));
		
		assert(subRate>=0 && subRate<=1) : "Substitution rate must be between 0 and 1, inclusive.  Invalid value: "+subRate;
		assert(indelRate>=0 && indelRate<=1) : "Indel rate must be between 0 and 1, inclusive.  Invalid value: "+indelRate;
		assert(errorRate>=0 && errorRate<=1) : "Total error rate must be between 0 and 1, inclusive.  Invalid value: "+errorRate;
		
		System.err.println(Tools.format("Target Identity:   \t%.2f%%", (1-errorRate2)*100));
		System.err.println(Tools.format("Substitution Rate: \t%.2f%%", subRate*100));
		System.err.println(Tools.format("Ins Rate:          \t%.2f%%", insRate*100));
		System.err.println(Tools.format("Del Rate:          \t%.2f%%", delRate*100));
		
		randy=Shared.threadLocalRandom(seed);
		
		{//Process parser fields
			Parser.processQuality();
			
			maxReads=parser.maxReads;
			in1=parser.in1;
			out1=parser.out1;
			overwrite=parser.overwrite;
			append=parser.append;
		}
		
		//Ensure there is an input file
		if(in1==null){throw new RuntimeException("Error - at least one input file is required.");}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1, outVcf)){
			outstream.println((out1==null)+", "+(outVcf==null)+", "+out1+", "+outVcf);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+", "+outVcf+"\n");
		}
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in1, out1, outVcf)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTA, null, true, overwrite, append, false);
		ffoutVcf=FileFormat.testOutput(outVcf, FileFormat.VCF, null, true, overwrite, append, false);
		ffin1=FileFormat.testInput(in1, FileFormat.FASTA, null, true, true);
	}
	
	private static int parsePloidy(String[] args){
		int ploidy=1;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(a.equals("ploidy")){
				ploidy=Integer.parseInt(b);
			}
		}
		return ploidy;
	}
	
	void process(Timer t){
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			cris.start();
		}

		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=4;
			
			if(cris.paired() && (in1==null || !in1.contains(".sam"))){
				outstream.println("Writing interleaved.");
			}

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, null, buff, null, false);
			ros.start();
		}else{ros=null;}

		{
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			if(reads!=null && !reads.isEmpty()){
				Read r=reads.get(0);
				assert((ffin1==null || ffin1.samOrBam()) || (r.mate!=null)==cris.paired());
			}
			ByteBuilder bb=new ByteBuilder();
			ArrayList<String> headers=(ffoutVcf==null ? null : new ArrayList<String>());
			ArrayList<SmallVar> vars=(ffoutVcf==null ? null : new ArrayList<SmallVar>());
			ArrayList<SmallVar> varsTemp=(ffoutVcf==null ? null : new ArrayList<SmallVar>());

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}
				
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);
					
					readsProcessed++;
					basesProcessed+=r1.length();
					
					processRead(r1, bb, varsTemp, headers);
					
					readsOut++;
					basesOut+=r1.length();
					if(vars!=null){vars.addAll(varsTemp);}
				}
				
				if(ros!=null){ros.add(reads, ln.id);}

				cris.returnList(ln);
				if(verbose){outstream.println("Returned a list.");}
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
			
			writeVars(vars, headers);
		}
		
		ReadWrite.closeStreams(cris, ros);
		if(verbose){outstream.println("Finished.");}
		
		
		{
			t.stop();
			
			//Calculate units per nanosecond
			double brnano=basesRetained/(double)(t.elapsed);
			
			//Add "k" and "m" for large numbers
			long mutationsAdded=subsAdded+insAdded+delsAdded+junctionsAdded;
			String brstring=(basesRetained<100000 ? ""+basesRetained : basesRetained<100000000 ? (basesRetained/1000)+"k" : (basesRetained/1000000)+"m");
			String mastring=(mutationsAdded<100000 ? ""+mutationsAdded : mutationsAdded<100000000 ? (mutationsAdded/1000)+"k" : (mutationsAdded/1000000)+"m");
			String rastring=(refAdded<100000 ? ""+refAdded : refAdded<100000000 ? (refAdded/1000)+"k" : (refAdded/1000000)+"m");
			String sastring=(subsAdded<100000 ? ""+subsAdded : subsAdded<100000000 ? (subsAdded/1000)+"k" : (subsAdded/1000000)+"m");
			String iastring=(insAdded<100000 ? ""+insAdded : insAdded<100000000 ? (insAdded/1000)+"k" : (insAdded/1000000)+"m");
			String dastring=(delsAdded<100000 ? ""+delsAdded : delsAdded<100000000 ? (delsAdded/1000)+"k" : (delsAdded/1000000)+"m");
			String ilastring=(insLenAdded<100000 ? ""+insLenAdded : insLenAdded<100000000 ? (insLenAdded/1000)+"k" : (insLenAdded/1000000)+"m");
			String dlastring=(delLenAdded<100000 ? ""+delLenAdded : delLenAdded<100000000 ? (delLenAdded/1000)+"k" : (delLenAdded/1000000)+"m");
			String nlastring=(netLengthAdded<100000 ? ""+netLengthAdded : netLengthAdded<100000000 ? (netLengthAdded/1000)+"k" : (netLengthAdded/1000000)+"m");
			String jastring=(junctionsAdded<100000 ? ""+junctionsAdded : junctionsAdded<100000000 ? (junctionsAdded/1000)+"k" : (junctionsAdded/1000000)+"m");
			
			//Format the strings so they have they are right-justified
			while(rastring.length()<8){rastring=" "+rastring;}
			while(ilastring.length()<8){ilastring=" "+ilastring;}
			while(dlastring.length()<8){dlastring=" "+dlastring;}
			while(nlastring.length()<8){nlastring=" "+nlastring;}
			
			while(brstring.length()<8){brstring=" "+brstring;}
			while(mastring.length()<8){mastring=" "+mastring;}
			while(sastring.length()<8){sastring=" "+sastring;}
			while(iastring.length()<8){iastring=" "+iastring;}
			while(dastring.length()<8){dastring=" "+dastring;}
			while(jastring.length()<8){jastring=" "+jastring;}

			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
			if(genomeFraction<1){outstream.println("Bases Retained:     "+brstring+" \t"+Tools.format("%.2fm bases/sec", brnano*1000));}
			outstream.println();
			outstream.println("Bases In:           "+Tools.padLeft(basesProcessed, 8));
			outstream.println("Bases Out:          "+Tools.padLeft(basesOut, 8));
			outstream.println("Mutations Added:    "+mastring+" \t"+Tools.format("%.2f%% Identity", 100f-mutationLengthAdded*100f/basesProcessed));
			outstream.println("Ref Added:          "+rastring);
			outstream.println("Subs Added:         "+sastring);
			outstream.println("Insertions Added:   "+iastring);
			outstream.println("Deletions Added:    "+dastring);
			outstream.println("Ins Length Added:   "+ilastring);
			outstream.println("Del Length Added:   "+dlastring);
			outstream.println("Net Indel Length:   "+nlastring);
			outstream.println("Junctions Added:    "+jastring);
		}
		
		t.stop();
	}
	
	public boolean isHomopolymerDel(byte[] bases, int pos, int len){
		final byte b=bases[pos];
		for(int i=1; i<len; i++){
			if(bases[pos+i]!=b){return false;}
		}
		if(pos>0 && bases[pos-1]==b){return true;}
		if(pos<bases.length-1 && bases[pos+1]==b){return true;}
		return false;
	}
	
	public boolean isHomopolymerIns(byte[] bases, int pos, byte b){
		if(b==bases[pos]){return true;}
		if(pos>0 && b==bases[pos-1]){return true;}
		return false;
	}
	
	public boolean isHomopolymerIns(byte[] bases, int pos, StringBuilder sb){
		byte b=(byte) sb.charAt(0);
		for(int i=1; i<sb.length(); i++) {
			if(sb.charAt(i)!=b){return false;}
		}
		return isHomopolymerIns(bases, pos, b);
	}
	
	public void processRead(Read r, ByteBuilder bb, ArrayList<SmallVar> vars, ArrayList<String> headers){
		
		if(r.aminoacid()) {
			processReadAmino(r, bb, vars, headers);
			return;
		}
		
		//Setup
		bb.clear();
		r.quality=null;
		if(headers!=null){headers.add("<ID="+r.id+",length="+r.length()+">");}
		if(vars!=null){vars.clear();}
		
		//Handle genomeFraction
		if(genomeFraction<1){
			final byte[] bases0=r.bases;
			int retain=(int)(bases0.length*(genomeFraction));
			if(retain<bases0.length){
				final int start=randy.nextInt(bases0.length);
				int i=0, j=start;
				for(; i<retain && j<bases0.length; i++, j++){
					bb.append(bases0[j]);
				}
				j=0;
				
				if(i<retain){
					junctionsAdded++;
					mutationLengthAdded++;
				} //Chimeric junction
				
				for(; i<retain; i++, j++){
					bb.append(bases0[j]);
				}
				r.bases=bb.toBytes();
				bb.clear();
			}
		}
		
		//Handle mutations
		final byte[] bases=r.bases;
		
		if(period>-1){
			int basesSinceMutation=0;
			char prevChar='N';
			for(int i=0; i<bases.length; i++){
				final byte b0=bases[i];
				byte b=b0;
				if(basesSinceMutation>=period && AminoAcid.isFullyDefined(b)){
					basesSinceMutation=0;
					float x=randy.nextFloat()*errorRate;
					if(x<subRate){
						b=AminoAcid.numberToBase[((AminoAcid.baseToNumber[b]+randy.nextInt(3)+1)&3)];
						bb.append(b);
						if(vars!=null){vars.add(new SmallVar(SUB, i, i+1, Character.toString((char)b), Character.toString((char)b0), prevChar, r.id, r.numericID));}
						subsAdded++;
						mutationLengthAdded++;
					}else if(randy.nextBoolean()){//del
						if(banHomopolymers && isHomopolymerDel(bases, i, 1)) {
							i--;
						}else{
							//do nothing
							if(vars!=null){vars.add(new SmallVar(DEL, i, i+1, "", Character.toString((char)b0), prevChar, r.id, r.numericID));}
							delsAdded++;
							delLenAdded++;
							mutationLengthAdded++;
							netLengthAdded--;
						}
					}else{//ins
						b=AminoAcid.numberToBase[randy.nextInt(4)];
						while(banHomopolymers && isHomopolymerIns(bases, i, b)) {
							b=AminoAcid.numberToBase[randy.nextInt(4)];
						}
						bb.append(b);
						if(vars!=null){vars.add(new SmallVar(INS, i, i, Character.toString((char)b), "", prevChar, r.id, r.numericID));}
						i--;
						insAdded++;
						insLenAdded++;
						mutationLengthAdded++;
						netLengthAdded++;
					}
				}else{
					basesSinceMutation++;
					refAdded++;
					bb.append(b);
				}
				prevChar=(char) b0;
			}
		}else{
			
			ConservationModel conservator=(sinewaves<1 ? null : 
				new ConservationModel(0.1f, sinewaves, randy));
			
			//All new code, rewritten from scratch, April 2025
			char prevChar='N';
			int lastIndel=-1;
			for(int i=0; i<bases.length;) {
			    final byte b0=bases[i];
			    final boolean defined=AminoAcid.isFullyDefined(b0);
			    final boolean spaceOK=(i-lastIndel>indelSpacing);
			    byte b=b0;
			    
			    if(conservator!=null && !conservator.shouldMutatePosition(i, randy)){
			        // Skip mutation, this is a conserved region
			        refAdded++;
			        bb.append(b);
			        i++;
			        continue;
			    }
			    
			    float x=randy.nextFloat();
			    boolean addSub=x<=subRate;
			    boolean addDel=x>subRate && x<=subRate+delRate && spaceOK;
			    boolean addIns=x>subRate+delRate && x<=errorRate && spaceOK;
			    boolean addRef=(i==0 || i==bases.length-1 || 
			    		x>errorRate || !defined || !(addDel || addIns || addSub));
			    boolean success=false;
			    if(bb.length()>0) {prevChar=(char) bb.get(bb.length()-1);}
//			    assert(false) : errorRate+", "+subRate+", "+delRate;
			    if(addRef) {
			        refAdded++;
			        bb.append(b);
			        success=true;
			        i++;
			    }else if(addSub) {
			    	subsAdded++;
			    	mutationLengthAdded++;
			    	b=AminoAcid.numberToBase[((AminoAcid.baseToNumber[b]+randy.nextInt(3)+1)&3)];
			    	bb.append(b);
			    	if(vars!=null){
			    		vars.add(new SmallVar(SUB, i, i+1, Character.toString((char)b), 
			    				Character.toString((char)b0), prevChar, r.id, r.numericID));
			    	}
			    	success=true;
			    	i++;
			    }else if(addIns) {
			    	int lim=Tools.min(maxIndel, bases.length-i-2), len=1;
					if(lim>=1){len=1+(Tools.min(randy.nextInt(lim), randy.nextInt(lim), randy.nextInt(lim)));}
//					if(insHomopolymer(r.bases, i, len, prevChar)) {
//			    		//Skip
//			    	}else 
			    	if(len>0) {
						insAdded++;
						insLenAdded+=len;
						mutationLengthAdded+=len;
						netLengthAdded+=len;
						ByteBuilder bbv=new ByteBuilder(len);
						for(int j=0; j<len; j++) {
					    	b=AminoAcid.numberToBase[randy.nextInt(4)];
					    	//b0 is the current character in the ref, which will be added after the insertion
					    	//prevChar is the last character in bb, which will come before the insertion
							while(banHomopolymers && ((j==0 && b==prevChar) || (j==len-1 && b==b0))) {
						    	b=AminoAcid.numberToBase[randy.nextInt(4)];
							}
					    	bb.append(b);
					    	bbv.append(b);
						}
						if(vars!=null){vars.add(new SmallVar(INS, i, i, bbv.toString(), "", prevChar, r.id, r.numericID));}
						success=true;
						lastIndel=i;
						//i++; //No i++ for insertions
					}
			    }else if(addDel) {
			    	int lim=Tools.min(maxIndel, bases.length-i-2), len=1;
			    	if(lim>=1){len=1+(Tools.min(randy.nextInt(lim), randy.nextInt(lim), randy.nextInt(lim)));}
			    	if(banHomopolymers && delHomopolymer(r.bases, i, len, prevChar)) {
			    		//Skip
			    	}else if(len>0) {
						delsAdded++;
						delLenAdded+=len;
						mutationLengthAdded+=len;
						netLengthAdded-=len;
						//do nothing
						if(vars!=null){vars.add(new SmallVar(DEL, i, i+len, "", 
								new String(bases, i, len), prevChar, r.id, r.numericID));}
						i=i+len;
						lastIndel=i;
						success=true;
			    	}
			    }
				
			    if(!success){
			    	// Problem encountered; advance
			    	refAdded++;
			    	bb.append(b);
					success=true;
			        i++;
			    }
			    
			    assert(success);
			    assert(bb.length()==refAdded+subsAdded+insLenAdded);
			    assert(i==refAdded+subsAdded+delLenAdded) : i+", "+refAdded+", "+subsAdded+", "+delLenAdded;
			}
		}
		
		if(padLeft>0) {
			ByteBuilder bb2=new ByteBuilder(bb.length+padLeft+padRight);
			for(int i=0; i<padLeft; i++) {bb2.append(AminoAcid.numberToBase[randy.nextInt(4)]);}
			bb2.append(bb);
			bb=bb2;
			if(vars!=null) {
				for(SmallVar v : vars) {
					v.start+=padLeft;
					v.stop+=padLeft;
				}
			}
		}
		for(int i=0; i<padRight; i++) {
			bb.append(AminoAcid.numberToBase[randy.nextInt(4)]);
		}
		
		condenseVars(vars);
		
		//Modify read
		r.bases=bb.toBytes();
		
		if(prefix!=null){
			r.id=prefix+r.numericID;
		}
		basesRetained+=r.bases.length;
	}
	
	private boolean insHomopolymer(byte[] bases, int i, int len, char prevChar) {
		return prevChar==bases[i];
	}
	
	private boolean delHomopolymer(byte[] bases, int i, int len, char prevChar) {
		byte a1=(byte)prevChar;
		byte a2=bases[i];
		byte b1=bases[i+len-1];
		byte b2=bases[i+len];
		
		if(a1==b2) {return true;}//It created a homopolymer (over-restrictive, but safe)
		if(a2==b1 && (a1==a2 || b1==b2)){return true;}//Partial deletion within a homopolymer
		return false;
	}

	public void processReadAmino(Read r, ByteBuilder bb, ArrayList<SmallVar> vars, ArrayList<String> headers){
		
		assert(r.aminoacid());
		
		//Setup
		bb.clear();
		r.quality=null;
		if(headers!=null){headers.add("<ID="+r.id+",length="+r.length()+">");}
		if(vars!=null){vars.clear();}
		
		//Handle genomeFraction
		if(genomeFraction<1){
			final byte[] bases0=r.bases;
			int retain=(int)(bases0.length*(genomeFraction));
			if(retain<bases0.length){
				final int start=randy.nextInt(bases0.length);
				int i=0, j=start;
				for(; i<retain && j<bases0.length; i++, j++){
					bb.append(bases0[j]);
				}
				j=0;
				
				if(i<retain){
					junctionsAdded++;
					mutationLengthAdded++;
				} //Chimeric junction
				
				for(; i<retain; i++, j++){
					bb.append(bases0[j]);
				}
				r.bases=bb.toBytes();
				bb.clear();
			}
		}
		
		//Handle mutations
		final byte[] bases=r.bases;
		
		if(period>-1){
			int basesSinceMutation=0;
			char prevChar='X';
			for(int i=0; i<bases.length; i++){
				final byte b0=bases[i];
				byte b=b0;
				if(basesSinceMutation>=period && AminoAcid.isFullyDefinedAA(b)){
					basesSinceMutation=0;
					float x=randy.nextFloat()*errorRate;
					if(x<subRate){
						subsAdded++;
						mutationLengthAdded++;
						b=AminoAcid.numberToAcid[((AminoAcid.acidToNumber[b]+randy.nextInt(20)+1)%21)];
						bb.append(b);
						if(vars!=null){vars.add(new SmallVar(SUB, i, i+1, Character.toString((char)b), Character.toString((char)b0), prevChar, r.id, r.numericID));}
					}else if(randy.nextBoolean()){//del
						delsAdded++;
						delLenAdded++;
						mutationLengthAdded++;
						netLengthAdded--;
						//do nothing
						if(vars!=null){vars.add(new SmallVar(DEL, i, i+1, "", Character.toString((char)b0), prevChar, r.id, r.numericID));}
					}else{//ins
						insAdded++;
						insLenAdded++;
						mutationLengthAdded++;
						netLengthAdded++;
						b=AminoAcid.numberToAcid[randy.nextInt(21)];
						bb.append(b);
						if(vars!=null){vars.add(new SmallVar(INS, i, i, Character.toString((char)b), "", prevChar, r.id, r.numericID));}
						i--;
					}
				}else{
					basesSinceMutation++;
					refAdded++;
					bb.append(b);
				}
				prevChar=(char) b0;
			}
		}else{
			char prevChar='N';
			for(int i=0; i<bases.length; i++){
				final byte b0=bases[i];
				byte b=b0;
				float x=randy.nextFloat();
				if(x<errorRate && AminoAcid.isFullyDefinedAA(b)){
					if(x<subRate){
						subsAdded++;
						mutationLengthAdded++;
						b=AminoAcid.numberToAcid[((AminoAcid.acidToNumber[b]+randy.nextInt(20)+1)%21)];
						bb.append(b);
						if(vars!=null){vars.add(new SmallVar(SUB, i, i+1, Character.toString((char)b), Character.toString((char)b0), prevChar, r.id, r.numericID));}
					}else if(randy.nextBoolean()){//del
						delsAdded++;
						delLenAdded++;
						mutationLengthAdded++;
						netLengthAdded--;
						//do nothing
						if(vars!=null){vars.add(new SmallVar(DEL, i, i+1, "", Character.toString((char)b0), prevChar, r.id, r.numericID));}
					}else{//ins
						insAdded++;
						insLenAdded++;
						mutationLengthAdded++;
						netLengthAdded++;
						b=AminoAcid.numberToAcid[randy.nextInt(21)];
						bb.append(b);
						if(vars!=null){vars.add(new SmallVar(INS, i, i, Character.toString((char)b), "", prevChar, r.id, r.numericID));}
						i--;
					}
				}else{
					refAdded++;
					bb.append(b);
				}
				prevChar=(char) b0;
			}
		}
		
		condenseVars(vars);
		
		//Modify read
		r.bases=bb.toBytes();
		
		if(prefix!=null){
			r.id=prefix+r.numericID;
		}
		basesRetained+=r.bases.length;
	}
	
	/*--------------------------------------------------------------*/
	
	private void condenseVars(ArrayList<SmallVar> vars){
		if(vars==null || vars.size()<2){return;}
		
		{//Pass 1: fuse indels into subs
			SmallVar current=null;
			for(int i=0; i<vars.size(); i++){
				SmallVar next=vars.get(i);
				if(next.type==SUB){
					current=null;
				}else if(current==null){
					current=next;
				}else if(current.stop==next.start && current.type!=next.type){
					if(current.type==DEL){
						assert(next.type==INS) : next.type;
						//Change the del to a sub
						current.type=SUB;
						current.alt=next.alt;
						current=null;
						vars.set(i, null);
					}else if(current.type==INS){
						assert(next.type==DEL) : next.type;
						//Change the ins to a sub
						current.type=SUB;
						current.ref=next.ref;
						current.stop=next.stop;
						current=null;
						vars.set(i, null);
					}else{
						assert(false) : current.type;
					}
				}else{
					current=next;
				}
			}
			Tools.condenseStrict(vars);
			if(vars.size()<2){return;}
		}
		
		{//Pass 2: lengthen indels
			SmallVar current=null;
			for(int i=0; i<vars.size(); i++){
				SmallVar next=vars.get(i);
				if(next.type==SUB){
					current=null;
					if(!next.valid()){vars.set(i, null);}
				}else if(current==null){
					current=next;
				}else if(current.stop==next.start && current.type==next.type){
					if(current.type==DEL){
						assert(next.type==DEL) : next.type;
						//Lengthen the deletion
						current.stop=next.stop;
						current.ref+=next.ref;
						vars.set(i, null);
					}else if(current.type==INS){
						assert(next.type==INS) : next.type;
						//Lengthen the insertion
						current.alt+=next.alt;
						vars.set(i, null);
					}else{
						assert(false) : current.type;
					}
				}else{
					current=next;
				}
			}
			Tools.condenseStrict(vars);
		}
	}
	
	void writeVars(ArrayList<SmallVar> vars, ArrayList<String> headers){
		if(ffoutVcf==null){return;}
		ByteStreamWriter bsw=new ByteStreamWriter(ffoutVcf);
		bsw.start();
		ByteBuilder bb=new ByteBuilder();
		bb.appendln("##fileformat=VCFv4.2");
		bb.appendln("##BBMapVersion="+Shared.BBTOOLS_VERSION_STRING);
		bb.appendln("##Program=MutateGenome");
		for(String s : headers){
			bb.append("##contig=").appendln(s);
		}
		bb.appendln("##FILTER=<ID=FAIL,Description=\"Fail\">");
		bb.appendln("##FILTER=<ID=PASS,Description=\"Pass\">");
		bb.appendln("##INFO=<ID=SN,Number=1,Type=Integer,Description=\"Scaffold Number\">");
		bb.appendln("##INFO=<ID=STA,Number=1,Type=Integer,Description=\"Start\">");
		bb.appendln("##INFO=<ID=STO,Number=1,Type=Integer,Description=\"Stop\">");
		bb.appendln("##INFO=<ID=TYP,Number=1,Type=String,Description=\"Type\">");
		bb.appendln("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
		bb.appendln("##FORMAT=<ID=SC,Number=1,Type=Float,Description=\"Score\">");
		bb.appendln("##FORMAT=<ID=PF,Number=1,Type=String,Description=\"Pass Filter\">");
		bb.appendln("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+(out1==null ? "sample" : ReadWrite.stripToCore(out1)));
		
		bsw.print(bb);
		bb.clear();
		
		for(SmallVar v : vars){
			v.toVcf(bb);
			bb.nl();
			if(bb.length()>=64000){
				bsw.print(bb);
				bb.clear();
			}
		}
		if(bb.length()>=0){
			bsw.print(bb);
			bb.clear();
		}
		errorState=bsw.poisonAndWait()|errorState;
	}
	
	/*--------------------------------------------------------------*/
	
	private static class SmallVar{
		
		SmallVar(int type_, int start_, int stop_, String alt_, String ref_, char prevChar_, String rname_, long scafNum_){
			type=type_;
			start=start_;
			stop=stop_;
			alt=alt_;
			ref=ref_;
			prevChar=prevChar_;
			rname=rname_;
			scafNum=scafNum_;
		}
		
		boolean valid(){
			return type!=SUB || !alt.equals(ref);
		}
		
		void toVcf(ByteBuilder bb){
			bb.append(rname).append('\t');
			if(type==SUB){
				bb.append(start+1).append('\t');
				bb.append('.').append('\t');
				bb.append(ref).append('\t');
				bb.append(alt).append('\t');
			}else if(type==DEL || type==INS){
				bb.append(start).append('\t');
				bb.append('.').append('\t');
				bb.append(prevChar).append(ref).append('\t');
				bb.append(prevChar).append(alt).append('\t');
			}else{assert(false);}
			bb.append("60.00").append('\t');
			bb.append("PASS").append('\t');
			bb.append("SN=").append(scafNum).append(';');
			bb.append("STA=").append(start).append(';');
			bb.append("STO=").append(stop).append(';');
			bb.append("TYP=").append(Var.typeArray[type]).append('\t');
			bb.append("GT:SC:PF").append('\t');
			bb.append(1).append(':');
			bb.append("60.00").append(':');
			bb.append("PASS");
		}
		
		int type;
		int start;
		int stop;
		String ref;
		String alt;
		final char prevChar;
		final String rname;
		final long scafNum;
		
	}
	
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String out1=null;
	private String outVcf=null;

	private String prefix=null;
	
	private final FileFormat ffin1;
	private final FileFormat ffout1;
	private final FileFormat ffoutVcf;
	
	/*--------------------------------------------------------------*/

	private long maxReads=-1;
//	private long mutationsAdded=0;

	private long mutationLengthAdded=0;
	private long netLengthAdded=0;
	private long subsAdded=0;
	private long refAdded=0;
	private long insAdded=0;
	private long delsAdded=0;
	private long insLenAdded=0;
	private long delLenAdded=0;
	private long junctionsAdded=0;

	int sinewaves=0;
	
	private int period=-1;
	
	private float genomeFraction=1;
	private long basesRetained;

	private long readsProcessed=0;
	private long basesProcessed=0;

	private long readsOut=0;
	private long basesOut=0;

	private float subRate=0;
	private float insRate=0;
	private float delRate=0;
	private int maxIndel=1;
	private int indelSpacing=3;
	private int padLeft=0, padRight=0;
	private boolean banHomopolymers=false;
	private final float errorRate;
	private final float errorRate2;
	
	private final Random randy;
	private long seed=-1;
	
	/*--------------------------------------------------------------*/
	
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;
	
	private static final int SUB=Var.SUB, INS=Var.INS, DEL=Var.DEL;
	
	private java.io.PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
