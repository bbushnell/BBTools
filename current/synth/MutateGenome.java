package synth;

import java.util.ArrayList;
import java.util.Arrays;
import shared.Random;

import bin.AdjustEntropy;
import bin.ConservationModel;
import clade.Clade;
import clade.CladeLoaderMF;
import clade.CladeObject;
import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
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

	/**
	 * Program entry point.
	 * Parses ploidy and delegates to MutateGenome2 for polyploid genomes.
	 * @param args Command-line arguments
	 */
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

	/**
	 * Constructs MutateGenome instance and parses command-line arguments.
	 * Configures mutation rates, file formats, and validation parameters.
	 * Sets up input/output streams and validates file accessibility.
	 * @param args Command-line arguments for configuration
	 */
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
			}else if(a.equalsIgnoreCase("preserveGC")){
				preserveGC=Parse.parseBoolean(b);
			}

			else if(a.equals("mod3")){
				mod3=Parse.parseBoolean(b);
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
			}

			else if(a.equals("rcomp")){//Doesn't do anything since it gets replaced.
				rcomp=Parse.parseBoolean(b);
			}else if(a.equals("counts")){//Doesn't do anything since it gets replaced.
				long[] counts=Parse.parseLongArray(b);
				kmerFreq=new float[counts.length];
				float mult=1f/Tools.sum(counts);
				for(int j=0; j<counts.length; k++) {kmerFreq[j]=counts[j]*mult;}
			}

			else if(a.equals("prefix")){
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
		if(k>0) {
			AdjustEntropy.load();//TODO: Annoying, should not be needed.
			Clade c=CladeLoaderMF.loadOneClade(in1, null, maxReads);
			kmerFreq=c.frequencies[k];
			rcomp=k>2;
		}
	}

	/**
	 * Parses ploidy value from command-line arguments.
	 * Used to determine whether to use MutateGenome or MutateGenome2.
	 * @param args Command-line arguments array
	 * @return Ploidy value (default 1)
	 */
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

	/**
	 * Main processing pipeline that applies mutations to input sequences.
	 * Reads sequences, applies configured mutations, and writes output files.
	 * Generates statistics and performance metrics upon completion.
	 * @param t Timer for tracking execution time
	 */
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
					final float gc=(preserveGC ? r1.gc() : 0.5f);

					readsProcessed++;
					basesProcessed+=r1.length();

					processRead(r1, bb, varsTemp, headers, gc);

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

	/**
	 * Checks if a deletion would create or extend a homopolymer run.
	 * Used to prevent homopolymer artifacts when banHomopolymers is enabled.
	 *
	 * @param bases The sequence bases array
	 * @param pos Starting position of potential deletion
	 * @param len Length of potential deletion
	 * @return true if deletion would create/extend homopolymer
	 */
	public boolean isHomopolymerDel(byte[] bases, int pos, int len){
		final byte b=bases[pos];
		for(int i=1; i<len; i++){
			if(bases[pos+i]!=b){return false;}
		}
		if(pos>0 && bases[pos-1]==b){return true;}
		if(pos<bases.length-1 && bases[pos+1]==b){return true;}
		return false;
	}

	/**
	 * Checks if an insertion would create or extend a homopolymer run.
	 *
	 * @param bases The sequence bases array
	 * @param pos Position for potential insertion
	 * @param b The base to insert
	 * @return true if insertion would create/extend homopolymer
	 */
	public boolean isHomopolymerIns(byte[] bases, int pos, byte b){
		if(b==bases[pos]){return true;}
		if(pos>0 && b==bases[pos-1]){return true;}
		return false;
	}

	/**
	 * Checks if a multi-base insertion would create homopolymer artifacts.
	 * Verifies the insertion sequence is homopolymeric before position check.
	 *
	 * @param bases The sequence bases array
	 * @param pos Position for potential insertion
	 * @param sb StringBuilder containing insertion sequence
	 * @return true if insertion would create/extend homopolymer
	 */
	public boolean isHomopolymerIns(byte[] bases, int pos, StringBuilder sb){
		byte b=(byte) sb.charAt(0);
		for(int i=1; i<sb.length(); i++) {
			if(sb.charAt(i)!=b){return false;}
		}
		return isHomopolymerIns(bases, pos, b);
	}

	public void processRead(Read r, ByteBuilder bb, ArrayList<SmallVar> vars, 
			ArrayList<String> headers, float gc){

		if(r.aminoacid()) {
			processReadAmino(r, bb, vars, headers);
			return;
		}

		long mutationLengthAdded=0;
		long netLengthAdded=0;
		long subsAdded=0;
		long refAdded=0;
		long insAdded=0;
		long delsAdded=0;
		long insLenAdded=0;
		long delLenAdded=0;
		long junctionsAdded=0;

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
						b=mutate(bases, i, gc, randy);
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
						b=AminoAcid.numberToBase[randy.nextInt()&3];
						while(banHomopolymers && isHomopolymerIns(bases, i, b)) {
							b=AminoAcid.numberToBase[randy.nextInt()&3];
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
					b=mutate(bases, i, gc, randy);
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
					if(mod3) {
						len=(len+2)%3;
						len=lim>len ? 0 : len;//TODO: Check if this is actually necessary
					}
					if(len>0) {
						insAdded++;
						insLenAdded+=len;
						mutationLengthAdded+=len;
						netLengthAdded+=len;
						ByteBuilder bbv=new ByteBuilder(len);
						for(int j=0; j<len; j++) {
							b=AminoAcid.numberToBase[randy.nextInt()&3];
							//b0 is the current character in the ref, which will be added after the insertion
							//prevChar is the last character in bb, which will come before the insertion
							while(banHomopolymers && ((j==0 && b==prevChar) || (j==len-1 && b==b0))) {
								b=AminoAcid.numberToBase[randy.nextInt()&3];
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
					if(mod3) {
						len=(len+2)%3;
						len=lim>len ? 0 : len;//TODO: Check if this is actually necessary
					}
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
				assert(bb.length()==refAdded+subsAdded+insLenAdded) : bb.length()+", "+refAdded+", "+
				subsAdded+", "+insLenAdded+", "+delLenAdded;//TODO: This is firing, claiming len 1.
				assert(i==refAdded+subsAdded+delLenAdded) : i+", "+refAdded+", "+subsAdded+", "+delLenAdded;
			}
		}

		this.mutationLengthAdded+=mutationLengthAdded;
		this.netLengthAdded+=netLengthAdded;
		this.subsAdded+=subsAdded;
		this.refAdded+=refAdded;
		this.insAdded+=insAdded;
		this.delsAdded+=delsAdded;
		this.insLenAdded+=insLenAdded;
		this.delLenAdded+=delLenAdded;
		this.junctionsAdded+=junctionsAdded;

		if(padLeft>0) {
			ByteBuilder bb2=new ByteBuilder(bb.length+padLeft+padRight);
			for(int i=0; i<padLeft; i++) {bb2.append(AminoAcid.numberToBase[randy.nextInt()&3]);}
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
			bb.append(AminoAcid.numberToBase[randy.nextInt()&3]);
		}

		condenseVars(vars);

		//Modify read
		r.bases=bb.toBytes();

		if(prefix!=null){
			r.id=prefix+r.numericID;
		}
		basesRetained+=r.bases.length;
	}
	
	
	public static byte mutate(final byte[] bases, final int pos, final float gc, Random randy) {
		final byte original=bases[pos];
		if(gc==0.5f) {
			int n=AminoAcid.baseToNumber[original];
			return AminoAcid.numberToBase[(1+randy.nextInt3()+n)&3];
		}
		
		// Determine if original is GC (0=A, 1=C, 2=G, 3=T)
		final boolean isGC = (original=='C' || original=='G'); // C or G
		
		// Decide whether to stay in same pool or switch
		final float stayProb = isGC ? gc : (1f - gc);
		final boolean stayInPool = randy.nextFloat() < stayProb;
		
		byte result;
		if(stayInPool) {
			// Stay in same GC/AT pool - pick the other base in that pool
			result=AminoAcid.baseToComplementExtended[original];
		} else {
			// Switch pools - pick randomly from the other pool
			boolean low=randy.nextBoolean();
			if(isGC) {
				result = low ? (byte)'A' : (byte)'T';
			} else {
				result = low ? (byte)'C' : (byte)'G';
			}
		}
		
		return result;
	}

	public static float[] calcPrefixProb(final byte[] bases, final int pos, final int k, final float[] freq, boolean rcomp) {
		float[] prob=new float[] {1, 1, 1, 1};
		final byte b0=bases[pos];
		final byte x0=AminoAcid.baseToNumber[b0];
		if(x0<0 || pos<k-1) {return prob;}
		int kmer=0;
		for(int i=0; i<k; i++) {
			final byte b=bases[pos-k+i+1];
			final byte x=AminoAcid.baseToNumber[b];
			if(x<0) {return prob;}
			kmer=((kmer<<2)|x);
		}
		kmer&=(~3);
		for(int i=0; i<4; i++) {
			int key=kmer|i;
			if(rcomp) {key=CladeObject.remapMatrix[k][key];}
			prob[i]=freq[key];
		}
		prob[x0]=0;
		return normalize(prob);
	}

	public static float[] normalize(float[] prob) {
		float sum=(float)Tools.sum(prob);
		if(sum==0) {
			Arrays.fill(prob, 1);
			sum=prob.length;
		}
		assert(sum!=0);
		return Tools.multiplyBy(prob, 1/sum);
	}

	/**
	 * Simplified homopolymer check for insertions.
	 *
	 * @param bases Sequence bases array
	 * @param i Position index
	 * @param len Insertion length (unused)
	 * @param prevChar Previous character in sequence
	 * @return true if insertion would extend homopolymer
	 */
	private boolean insHomopolymer(byte[] bases, int i, int len, char prevChar) {
		return prevChar==bases[i];
	}

	/**
	 * Checks if deletion would create homopolymer artifacts.
	 * Examines bases before, within, and after deletion region.
	 *
	 * @param bases Sequence bases array
	 * @param i Starting position of deletion
	 * @param len Length of deletion
	 * @param prevChar Previous character in output sequence
	 * @return true if deletion would create homopolymer artifacts
	 */
	private boolean delHomopolymer(byte[] bases, int i, int len, char prevChar) {
		byte a1=(byte)prevChar;
		byte a2=bases[i];
		byte b1=bases[i+len-1];
		byte b2=bases[i+len];

		if(a1==b2) {return true;}//It created a homopolymer (over-restrictive, but safe)
		if(a2==b1 && (a1==a2 || b1==b2)){return true;}//Partial deletion within a homopolymer
		return false;
	}

	/**
	 * Processes amino acid sequences with mutations specific to protein sequences.
	 * Uses 21-amino acid alphabet and appropriate mutation logic.
	 *
	 * @param r The amino acid read to mutate
	 * @param bb ByteBuilder for sequence construction
	 * @param vars ArrayList to collect variant information
	 * @param headers ArrayList to collect VCF headers
	 */
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

	/**
	 * Condenses adjacent variants into more complex variants.
	 * Two-pass algorithm: first fuses indels into substitutions,
	 * then lengthens consecutive indels of same type.
	 * @param vars ArrayList of variants to condense
	 */
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

	/**
	 * Writes variants to VCF format output file.
	 * Includes standard VCF headers and variant records.
	 * @param vars ArrayList of variants to write
	 * @param headers ArrayList of contig headers
	 */
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

		/**
		 * Constructs a variant record with position and sequence information.
		 *
		 * @param type_ Variant type (SUB, INS, or DEL)
		 * @param start_ Start position in sequence
		 * @param stop_ Stop position in sequence
		 * @param alt_ Alternative allele sequence
		 * @param ref_ Reference allele sequence
		 * @param prevChar_ Character preceding variant position
		 * @param rname_ Reference sequence name
		 * @param scafNum_ Scaffold numeric identifier
		 */
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

		/** Validates that substitution variants have different ref and alt alleles.
		 * @return true if variant is valid, false if ref equals alt for substitution */
		boolean valid(){
			return type!=SUB || !alt.equals(ref);
		}

		/**
		 * Formats variant as VCF record line.
		 * Handles coordinate adjustments for different variant types.
		 * @param bb ByteBuilder to append VCF format line
		 */
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

		/** Variant type (SUB, INS, or DEL) */
		int type;
		/** Start position of variant in sequence */
		int start;
		/** Stop position of variant in sequence */
		int stop;
		/** Reference allele sequence */
		String ref;
		/** Alternative allele sequence */
		String alt;
		/** Character preceding variant position */
		final char prevChar;
		/** Reference sequence name */
		final String rname;
		/** Numeric scaffold identifier */
		final long scafNum;

	}

	/*--------------------------------------------------------------*/

	/** Input file path */
	private String in1=null;
	/** Output FASTA file path */
	private String out1=null;
	/** Output VCF file path */
	private String outVcf=null;

	/** Prefix to apply to sequence names */
	private String prefix=null;

	/** Input file format specification */
	private final FileFormat ffin1;
	/** Output FASTA file format specification */
	private final FileFormat ffout1;
	/** Output VCF file format specification */
	private final FileFormat ffoutVcf;

	/*--------------------------------------------------------------*/

	/** Maximum number of reads to process (-1 for unlimited) */
	private long maxReads=-1;

	/** Total length of mutations added to sequences */
	private long mutationLengthAdded=0;
	/** Net length change from insertions minus deletions */
	private long netLengthAdded=0;
	/** Count of substitution mutations added */
	private long subsAdded=0;
	/** Count of reference bases retained without mutation */
	private long refAdded=0;
	/** Count of insertion mutations added */
	private long insAdded=0;
	/** Count of deletion mutations added */
	private long delsAdded=0;
	/** Total length of all insertions added */
	private long insLenAdded=0;
	/** Total length of all deletions added */
	private long delLenAdded=0;
	/** Count of chimeric junctions from genome fraction sampling */
	private long junctionsAdded=0;

	/** Number of sine waves for conservation model (0=disabled) */
	int sinewaves=0;
	boolean preserveGC=true;

	/** Period for regular mutation pattern (-1 for random) */
	private int period=-1;

	/** Fraction of genome to retain (1.0 = complete genome) */
	private float genomeFraction=1;
	/** Total bases retained in output sequences */
	private long basesRetained;

	/** Count of input sequences processed */
	private long readsProcessed=0;
	/** Total bases processed from input sequences */
	private long basesProcessed=0;

	/** Count of output sequences written */
	private long readsOut=0;
	/** Total bases written to output sequences */
	private long basesOut=0;

	/** Rate of substitution mutations (0.0-1.0) */
	private float subRate=0;
	/** Rate of insertion mutations (0.0-1.0) */
	private float insRate=0;
	/** Rate of deletion mutations (0.0-1.0) */
	private float delRate=0;
	/** Maximum length for insertion or deletion mutations */
	private int maxIndel=1;
	/** Minimum spacing between adjacent indel mutations */
	private int indelSpacing=3;
	/** Number of random bases to add at sequence start */
	private int padLeft=0, padRight=0;
	/** Whether to prevent mutations that create homopolymer runs */
	private boolean banHomopolymers=false;
	/** Total error rate combining all mutation types */
	private final float errorRate;
	/** Adjusted error rate accounting for indel length distribution */
	private final float errorRate2;

	private boolean mod3=false;
	private float[] kmerFreq;
	private int k=-1;
	private boolean rcomp=false;

	/** Random number generator for mutation decisions */
	private final Random randy;
	/** Random seed for reproducible mutation patterns */
	private long seed=-1;

	/*--------------------------------------------------------------*/

	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=true;
	/** Append to existing output files */
	private boolean append=false;

	/** Constant for substitution variant type */
	private static final int SUB=Var.SUB, INS=Var.INS, DEL=Var.DEL;

	/** Output stream for status and error messages */
	private java.io.PrintStream outstream=System.err;
	/** Enable verbose logging output */
	public static boolean verbose=false;

}
