package synth;

import java.util.ArrayList;
import java.util.Arrays;
import shared.Random;

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
 * Polyploid synthetic genome mutator used by {@link MutateGenome} when {@code ploidy!=1}.
 * Emits one mutated FASTA record per input record for haploid mode, or one record per
 * haplotype for polyploid mode, and can optionally record generated variants as VCF.
 * Supports both nucleotide and amino-acid inputs; amino records use the alternate mutation
 * path selected by {@link Read#aminoacid()}.
 *
 * @author Brian Bushnell
 * @date June 1, 2016
 */
public class MutateGenome2 {

	public static void main(String[] args){
		Timer t=new Timer();
		MutateGenome2 x=new MutateGenome2(args);
		x.process(t);

		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}

	/**
	 * Parses parameters, reports the requested aggregate mutation rates, initializes the
	 * random stream and binds FASTA/VCF input and output formats.
	 */
	public MutateGenome2(String[] args){

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
				indelRate=Float.parseFloat(b);
				if(indelRate>1){indelRate/=100;}
			}else if(a.equals("maxindel")){
				maxIndel=Parse.parseIntKMG(b);
			}else if(a.equals("indelspacing")){
				indelSpacing=Parse.parseIntKMG(b);
			}else if(a.equals("seed")){
				seed=Long.parseLong(b);
			}else if(a.equals("ploidy")){
				ploidy=Integer.parseInt(b);
			}else if(a.equals("nohomopolymers") || a.equals("nohomo") || a.equals("banhomopolymers") || a.equals("banhomo")){
				banHomopolymers=Parse.parseBoolean(b);
			}else if(a.equalsIgnoreCase("hetrate")){
				hetRate=Float.parseFloat(b);
			}else if(a.equals("prefix")){
				prefix=b;
			}else if(a.equals("vcf") || a.equals("outvcf") || a.equals("vcfout") || a.equals("vars") || a.equals("varsout") || a.equals("outvars")){
				outVcf=b;
			}else if(a.equals("id") || a.equals("identity")){
				float x=Float.parseFloat(b);
				if(x>1){x=x/100;}
				x=1-x;
				indelRate=x*0.01f;
				subRate=x*0.99f;
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
		errorRate=subRate+indelRate;
		errorRate2=subRate+(indelRate*0.16666f*(maxIndel+1));

		assert(subRate>=0 && subRate<=1) : "Substitution rate must be between 0 and 1, inclusive.  Invalid value: "+subRate;
		assert(indelRate>=0 && indelRate<=1) : "Indel rate must be between 0 and 1, inclusive.  Invalid value: "+indelRate;
		assert(errorRate>=0 && errorRate<=1) : "Total error rate must be between 0 and 1, inclusive.  Invalid value: "+errorRate;

		System.err.println(Tools.format("Target Identity:   \t%.2f%%", (1-errorRate2)*100));
		System.err.println(Tools.format("Substitution Rate: \t%.2f%%", subRate*100));
		System.err.println(Tools.format("Indel Rate:        \t%.2f%%", indelRate*100));

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

	/**
	 * Streams input records, expands each into one or more mutated haplotypes, writes the
	 * generated FASTA records, then emits the accumulated VCF variants.
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
			ByteBuilder[] bba=new ByteBuilder[ploidy];
			for(int i=0; i<bba.length; i++){bba[i]=new ByteBuilder();}
			ArrayList<String> headers=(ffoutVcf==null ? null : new ArrayList<String>());
			ArrayList<SmallVar> vars=(ffoutVcf==null ? null : new ArrayList<SmallVar>());
			ArrayList<SmallVar> varsTemp=(ffoutVcf==null ? null : new ArrayList<SmallVar>());

			while(ln!=null && reads!=null && reads.size()>0){//ln!=null prevents a compiler potential null access warning
				if(verbose){outstream.println("Fetched "+reads.size()+" reads.");}

				ArrayList<Read> mutants=new ArrayList<Read>(ploidy*reads.size());
				for(int idx=0; idx<reads.size(); idx++){
					final Read r1=reads.get(idx);

					readsProcessed++;
					basesProcessed+=r1.length();

					ArrayList<Read> mutant=processRead(r1, bba, varsTemp, headers);
					mutants.addAll(mutant);
					if(vars!=null){vars.addAll(varsTemp);}
				}

				if(ros!=null){ros.add(mutants, ln.id);}

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

		errorState|=ReadWrite.closeStreams(cris, ros);
		if(verbose){outstream.println("Finished.");}
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}


		{
			t.stop();

			//Calculate units per nanosecond
			double brnano=basesRetained/(double)(t.elapsed);

			//Add "k" and "m" for large numbers
			long mutationsAdded=subsAdded+insAdded+delsAdded+junctionsAdded;
			String brstring=(basesRetained<100000 ? ""+basesRetained : basesRetained<100000000 ? (basesRetained/1000)+"k" : (basesRetained/1000000)+"m");
			String mastring=(mutationsAdded<100000 ? ""+mutationsAdded : mutationsAdded<100000000 ? (mutationsAdded/1000)+"k" : (mutationsAdded/1000000)+"m");
			String sastring=(subsAdded<100000 ? ""+subsAdded : subsAdded<100000000 ? (subsAdded/1000)+"k" : (subsAdded/1000000)+"m");
			String iastring=(insAdded<100000 ? ""+insAdded : insAdded<100000000 ? (insAdded/1000)+"k" : (insAdded/1000000)+"m");
			String dastring=(delsAdded<100000 ? ""+delsAdded : delsAdded<100000000 ? (delsAdded/1000)+"k" : (delsAdded/1000000)+"m");
			String jastring=(junctionsAdded<100000 ? ""+junctionsAdded : junctionsAdded<100000000 ? (junctionsAdded/1000)+"k" : (junctionsAdded/1000000)+"m");

			//Format the strings so they have they are right-justified
			while(brstring.length()<8){brstring=" "+brstring;}
			while(mastring.length()<8){mastring=" "+mastring;}
			while(sastring.length()<8){sastring=" "+sastring;}
			while(iastring.length()<8){iastring=" "+iastring;}
			while(dastring.length()<8){dastring=" "+dastring;}
			while(jastring.length()<8){jastring=" "+jastring;}

			outstream.println(Tools.timeReadsBasesProcessed(t, readsProcessed, basesProcessed, 8));
			if(genomeFraction<1){outstream.println("Bases Retained:     "+brstring+" \t"+Tools.format("%.2fm bases/sec", brnano*1000));}
			outstream.println("Mutations Added:    "+mastring+" \t"+Tools.format("%.2f%% Identity", 100f-mutationLengthAdded*100f/basesProcessed));
			outstream.println("Subs Added:         "+sastring);
			outstream.println("Insertions Added:   "+iastring);
			outstream.println("Deletions Added:    "+dastring);
			outstream.println("Junctions Added:    "+jastring);
		}

		t.stop();
	}

	/**
	 * Creates the per-haplotype presence mask for the next variant.
	 * A value of 1 means the mutation is present on that haplotype; 0 keeps the reference.
	 * Heterozygous calls require at least one mutated and one reference haplotype.
	 */
	public int[] makePresentArray(){
		int[] present=new int[ploidy];
		if(ploidy>1 && randy.nextFloat()<hetRate){
			int sum=0;
			while(sum<1 || sum>=ploidy){
				sum=0;
				Arrays.fill(present, 0);
				for(int i=0; i<ploidy; i++){
					int x=randy.nextInt2();
					sum+=x;
					present[i]=x;
				}
			}
		}else{
			Arrays.fill(present, 1);
		}
		return present;
	}

	public boolean isHomopolymerDel(byte[] bases, int pos, int len){
		final byte b=bases[pos];
		for(int i=1; i<len; i++){
			if(bases[pos+i]!=b){return false;}
		}
		if(pos>0 && bases[pos-1]==b){return true;}
		if(pos+len<bases.length && bases[pos+len]==b){return true;}
		return false;
	}

	/** Returns true if inserting {@code b} at {@code pos} would create a homopolymer. */
	public boolean isHomopolymerIns(byte[] bases, int pos, byte b){
		if(b==bases[pos]){return true;}
		if(pos>0 && b==bases[pos-1]){return true;}
		return false;
	}

	/**
	 * Returns true if the proposed insertion is a single repeated nucleotide and that base
	 * would create a homopolymer at {@code pos}.
	 */
	public boolean isHomopolymerIns(byte[] bases, int pos, StringBuilder sb){
		byte b=(byte) sb.charAt(0);
		for(int i=1; i<sb.length(); i++) {
			if(sb.charAt(i)!=b){return false;}
		}
		return isHomopolymerIns(bases, pos, b);
	}

	/**
	 * Mutates one nucleotide record, optionally wraps it to a random retained fraction, and
	 * returns one output read per haplotype. {@code vars} is cleared and filled with this
	 * record's variants; the VCF contig header is appended when VCF output is enabled.
	 */
	public ArrayList<Read> processRead(Read r, ByteBuilder[] bba, ArrayList<SmallVar> vars, ArrayList<String> headers){

		if(r.aminoacid()) {
			return processReadAmino(r, bba, vars, headers);
		}

		//Setup
		for(ByteBuilder bb : bba){bb.clear();}
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
					bba[0].append(bases0[j]);
				}
				j=0;

				if(i<retain){
					junctionsAdded++;
					mutationLengthAdded++;
				} //Chimeric junction

				for(; i<retain; i++, j++){
					bba[0].append(bases0[j]);
				}
				r.bases=bba[0].toBytes();
				bba[0].clear();
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
					int[] present=makePresentArray();
					if(x<subRate){
						b=AminoAcid.numberToBase[((AminoAcid.baseToNumber[b]+randy.nextInt3()+1)&3)];

						for(int haplo=0; haplo<ploidy; haplo++){
							if(present[haplo]==1){
								bba[haplo].append(b);
							}else{
								bba[haplo].append(b0);
							}
						}
						if(vars!=null){vars.add(new SmallVar(SUB, i, i+1, Character.toString((char)b), Character.toString((char)b0), prevChar, r.id, r.numericID, present));}
						subsAdded++;
						mutationLengthAdded++;
					}else if(randy.nextBoolean()){//del
						if(banHomopolymers && isHomopolymerDel(bases, i, 1)) {
							i--;
						}else{
							for(int haplo=0; haplo<ploidy; haplo++){
								if(present[haplo]==1){
									//do nothing
								}else{
									bba[haplo].append(b0);
								}
							}
							if(vars!=null){vars.add(new SmallVar(DEL, i, i+1, "", Character.toString((char)b0), prevChar, r.id, r.numericID, present));}
							delsAdded++;
							mutationLengthAdded++;
						}
					}else{//ins
						b=AminoAcid.numberToBase[randy.nextInt()&3];
						while(banHomopolymers && isHomopolymerIns(bases, i, b)) {
							b=AminoAcid.numberToBase[randy.nextInt()&3];
						}
						for(int haplo=0; haplo<ploidy; haplo++){
							if(present[haplo]==1){
								bba[haplo].append(b);
							}else{
								//do nothing
							}
						}
						if(vars!=null){vars.add(new SmallVar(INS, i, i, Character.toString((char)b), "", prevChar, r.id, r.numericID, present));}
						i--;
						insAdded++;
						mutationLengthAdded++;
					}
				}else{
					basesSinceMutation++;
					for(ByteBuilder bb : bba){
						bb.append(b);
					}
				}
				prevChar=(char) b0;
			}
		}else{
			char prevChar='N';
			int lastIndel=-1;
			for(int i=0; i<bases.length; i++){
				final byte b0=bases[i];
				byte b=b0;
				float x=randy.nextFloat();
				if(x<errorRate && AminoAcid.isFullyDefined(b)){
					int[] present=makePresentArray();
					//					System.err.println(x+", "+errorRate+", "+subRate);
					if(x<subRate){
						subsAdded++;
						mutationLengthAdded++;
						b=AminoAcid.numberToBase[((AminoAcid.baseToNumber[b]+randy.nextInt3()+1)&3)];
						for(int haplo=0; haplo<ploidy; haplo++){
							if(present[haplo]==1){
								bba[haplo].append(b);
							}else{
								bba[haplo].append(b0);
							}
						}
						if(vars!=null){vars.add(new SmallVar(SUB, i, i+1, Character.toString((char)b), Character.toString((char)b0), prevChar, r.id, r.numericID, present));}
					}else{
						boolean attemptedIndel=false;
						if(i-lastIndel>indelSpacing){
							int lim=Tools.min(maxIndel, bases.length-i-1);
							if(lim>=1){
								attemptedIndel=true;
								int len=1+(Tools.min(randy.nextInt(lim), randy.nextInt(lim), randy.nextInt(lim)));

								if(randy.nextBoolean()){//del
									if(banHomopolymers && isHomopolymerDel(bases, i, len)) {
										i--;
									}else{
										delsAdded++;
										mutationLengthAdded+=len;
										byte[] ref=new byte[len];
										for(int rpos=0; rpos<len; rpos++){ref[rpos]=bases[i+rpos];}
										//do nothing
										for(int haplo=0; haplo<ploidy; haplo++){
											if(present[haplo]==1){
												//do nothing
											}else{
												bba[haplo].append(ref);
											}
										}
										if(vars!=null){vars.add(new SmallVar(DEL, i, i+len, "", new String(bases, i, len), prevChar, r.id, r.numericID, present));}
										i=i+len-1;
									}
								}else{//ins
									insAdded++;
									mutationLengthAdded+=len;
									if(len==1){
										b=AminoAcid.numberToBase[randy.nextInt()&3];
										while(banHomopolymers && isHomopolymerIns(bases, i, b)) {
											b=AminoAcid.numberToBase[randy.nextInt()&3];
										}
										for(int haplo=0; haplo<ploidy; haplo++){
											if(present[haplo]==1){
												bba[haplo].append(b);
											}else{
												//do nothing
											}
										}
										if(vars!=null){vars.add(new SmallVar(INS, i, i, Character.toString((char)b), "", prevChar, r.id, r.numericID, present));}
									}else{
										StringBuilder sb=new StringBuilder();
										while(sb.length()==0 || (banHomopolymers && isHomopolymerIns(bases, i, sb))) {
											sb.setLength(0);
											while(sb.length()<len){
												b=AminoAcid.numberToBase[randy.nextInt()&3];
												sb.append((char)b);
											}
										}
										for(int haplo=0; haplo<ploidy; haplo++){
											if(present[haplo]==1){
												bba[haplo].append(sb);
											}else{
												//do nothing
											}
										}
										if(vars!=null){vars.add(new SmallVar(INS, i, i, sb.toString(), "", prevChar, r.id, r.numericID, present));}
									}
									i--;
								}
								lastIndel=i;
							}
						}
						if(!attemptedIndel){
							for(ByteBuilder bb : bba){
								bb.append(b);
							}
						}
					}
				}else{
					for(ByteBuilder bb : bba){
						bb.append(b);
					}
				}
				prevChar=(char) b0;
			}
		}

		condenseVars(vars);

		//Modify read
		r.bases=bba[0].toBytes();

		if(prefix!=null){
			r.id=prefix+r.numericID;
		}
		basesRetained+=r.bases.length;

		ArrayList<Read> ret=new ArrayList<Read>();
		if(ploidy==1){ret.add(r);}
		else{
			for(int i=0; i<ploidy; i++){
				ret.add(new Read(bba[i].toBytes(), null, r.id+"_haplo_"+i, r.numericID));
			}
		}
		return ret;
	}

	/**
	 * Amino-acid analogue of {@link #processRead}; uses amino-acid alphabets and skips the
	 * nucleotide homopolymer guard.
	 */
	public ArrayList<Read> processReadAmino(Read r, ByteBuilder[] bba, ArrayList<SmallVar> vars, ArrayList<String> headers){

		assert(r.aminoacid());

		//Setup
		for(ByteBuilder bb : bba){bb.clear();}
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
					bba[0].append(bases0[j]);
				}
				j=0;

				if(i<retain){
					junctionsAdded++;
					mutationLengthAdded++;
				} //Chimeric junction

				for(; i<retain; i++, j++){
					bba[0].append(bases0[j]);
				}
				r.bases=bba[0].toBytes();
				bba[0].clear();
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
					int[] present=makePresentArray();
					if(x<subRate){
						b=AminoAcid.numberToAcid[((AminoAcid.acidToNumber[b]+randy.nextInt(20)+1)%21)];

						for(int haplo=0; haplo<ploidy; haplo++){
							if(present[haplo]==1){
								bba[haplo].append(b);
							}else{
								bba[haplo].append(b0);
							}
						}
						if(vars!=null){vars.add(new SmallVar(SUB, i, i+1, Character.toString((char)b), Character.toString((char)b0), prevChar, r.id, r.numericID, present));}
						subsAdded++;
						mutationLengthAdded++;
					}else if(randy.nextBoolean()){//del
						for(int haplo=0; haplo<ploidy; haplo++){
							if(present[haplo]==1){
								//do nothing
							}else{
								bba[haplo].append(b0);
							}
						}
						if(vars!=null){vars.add(new SmallVar(DEL, i, i+1, "", Character.toString((char)b0), prevChar, r.id, r.numericID, present));}
						delsAdded++;
						mutationLengthAdded++;
					}else{//ins
						b=AminoAcid.numberToAcid[randy.nextInt(21)];
						for(int haplo=0; haplo<ploidy; haplo++){
							if(present[haplo]==1){
								bba[haplo].append(b);
							}else{
								//do nothing
							}
						}
						if(vars!=null){vars.add(new SmallVar(INS, i, i, Character.toString((char)b), "", prevChar, r.id, r.numericID, present));}
						i--;
						insAdded++;
						mutationLengthAdded++;
					}
				}else{
					basesSinceMutation++;
					for(ByteBuilder bb : bba){
						bb.append(b);
					}
				}
				prevChar=(char) b0;
			}
		}else{
			char prevChar='N';
			int lastIndel=-1;
			for(int i=0; i<bases.length; i++){
				final byte b0=bases[i];
				byte b=b0;
				float x=randy.nextFloat();
				if(x<errorRate && AminoAcid.isFullyDefinedAA(b)){
					int[] present=makePresentArray();
					//					System.err.println(x+", "+errorRate+", "+subRate);
					if(x<subRate){
						subsAdded++;
						mutationLengthAdded++;
						b=AminoAcid.numberToAcid[((AminoAcid.acidToNumber[b]+randy.nextInt(20)+1)%21)];
						for(int haplo=0; haplo<ploidy; haplo++){
							if(present[haplo]==1){
								bba[haplo].append(b);
							}else{
								bba[haplo].append(b0);
							}
						}
						if(vars!=null){vars.add(new SmallVar(SUB, i, i+1, Character.toString((char)b), Character.toString((char)b0), prevChar, r.id, r.numericID, present));}
					}else{
						boolean attemptedIndel=false;
						if(i-lastIndel>indelSpacing){
							int lim=Tools.min(maxIndel, bases.length-i-1);
							if(lim>=1){
								attemptedIndel=true;
								int len=1+(Tools.min(randy.nextInt(lim), randy.nextInt(lim), randy.nextInt(lim)));

								if(randy.nextBoolean()){//del
									delsAdded++;
									mutationLengthAdded+=len;
									byte[] ref=new byte[len];
									for(int rpos=0; rpos<len; rpos++){ref[rpos]=bases[i+rpos];}
									//do nothing
									for(int haplo=0; haplo<ploidy; haplo++){
										if(present[haplo]==1){
											//do nothing
										}else{
											bba[haplo].append(ref);
										}
									}
									if(vars!=null){vars.add(new SmallVar(DEL, i, i+len, "", new String(bases, i, len), prevChar, r.id, r.numericID, present));}
									i=i+len-1;
								}else{//ins
									insAdded++;
									mutationLengthAdded+=len;
									if(len==1){
										b=AminoAcid.numberToAcid[randy.nextInt(21)];
										for(int haplo=0; haplo<ploidy; haplo++){
											if(present[haplo]==1){
												bba[haplo].append(b);
											}else{
												//do nothing
											}
										}
										if(vars!=null){vars.add(new SmallVar(INS, i, i, Character.toString((char)b), "", prevChar, r.id, r.numericID, present));}
									}else{
										StringBuilder sb=new StringBuilder();
										while(sb.length()<len){
											b=AminoAcid.numberToAcid[randy.nextInt(21)];
											sb.append((char)b);
										}
										for(int haplo=0; haplo<ploidy; haplo++){
											if(present[haplo]==1){
												bba[haplo].append(sb);
											}else{
												//do nothing
											}
										}
										if(vars!=null){vars.add(new SmallVar(INS, i, i, sb.toString(), "", prevChar, r.id, r.numericID, present));}
									}
									i--;
								}
								lastIndel=i;
							}
						}
						if(!attemptedIndel){
							for(int haplo=0; haplo<ploidy; haplo++){
								if(present[haplo]==1){
									bba[haplo].append(b);
								}else{
									bba[haplo].append(b0);
								}
							}
						}
					}
				}else{
					for(ByteBuilder bb : bba){
						bb.append(b);
					}
				}
				prevChar=(char) b0;
			}
		}

		condenseVars(vars);

		//Modify read
		r.bases=bba[0].toBytes();

		if(prefix!=null){
			r.id=prefix+r.numericID;
		}
		basesRetained+=r.bases.length;

		ArrayList<Read> ret=new ArrayList<Read>();
		if(ploidy==1){ret.add(r);}
		else{
			for(int i=0; i<ploidy; i++){
				ret.add(new Read(bba[i].toBytes(), null, r.id+"_haplo_"+i, r.numericID));
			}
		}
		return ret;
	}

	/*--------------------------------------------------------------*/

	/**
	 * Combines adjacent same-phase indels where the simplified VCF can represent them as a
	 * substitution or as one longer same-type indel.
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
				}else if(current.stop==next.start && current.type!=next.type && current.sharesPhase(next)){
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
				}else if(current.stop==next.start && current.type==next.type && current.sharesPhase(next)){
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

	/** Writes generated variants after input records are processed and FASTA batches enqueued. */
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

	/** Minimal zero-based variant record plus the haplotypes on which it is present. */
	private static class SmallVar{

		SmallVar(int type_, int start_, int stop_, String alt_, String ref_, char prevChar_, String rname_, long scafNum_, int[] present_){
			type=type_;
			start=start_;
			stop=stop_;
			alt=alt_;
			ref=ref_;
			prevChar=prevChar_;
			rname=rname_;
			scafNum=scafNum_;
			present=present_;
		}

		/** Returns false only for condensed substitutions that did not change bases. */
		boolean valid(){
			return type!=SUB || !alt.equals(ref);
		}

		/** Appends one VCF data row, converting to VCF's coordinate and indel anchoring rules. */
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
			bb.append(present[0]);
			for(int i=1; i<present.length; i++){bb.append('|').append(present[i]);}
			bb.append(':');
			bb.append("60.00").append(':');
			bb.append("PASS");
		}

		/** Returns true when both variants are present on exactly the same haplotypes. */
		boolean sharesPhase(SmallVar sv){
			for(int i=0; i<present.length; i++){
				if(present[i]!=sv.present[i]){return false;}
			}
			return true;
		}

		/** Variant type from {@link Var}. */
		int type;
		/** Zero-based inclusive start for insertions, zero-based start for deletions/subs. */
		final int start;
		/** Zero-based exclusive stop; equal to start for insertions. */
		int stop;
		/** Reference allele without the previous-base anchor. */
		String ref;
		/** Alternate allele without the previous-base anchor. */
		String alt;
		/** Base before the event, used as the VCF anchor for indels. */
		final char prevChar;
		/** Input record ID used as the VCF CHROM field. */
		final String rname;
		/** Input record numeric ID emitted as the SN info field. */
		final long scafNum;

		/** One element per haplotype: 1 for mutated, 0 for reference. */
		final int[] present;

	}

	/*--------------------------------------------------------------*/

	/** Input FASTA path. */
	private String in1=null;
	/** Mutated FASTA output path. */
	private String out1=null;
	/** Optional VCF output path. */
	private String outVcf=null;

	/** Optional output ID prefix; replaces the original ID with prefix+numericID. */
	private String prefix=null;

	/** Validated input FASTA format. */
	private final FileFormat ffin1;
	/** Validated output FASTA format. */
	private final FileFormat ffout1;
	/** Validated optional VCF format. */
	private final FileFormat ffoutVcf;

	/*--------------------------------------------------------------*/

	private long maxReads=-1;
//	private long mutationsAdded=0;

	/** Number of haplotypes to emit per input record. */
	private int ploidy=1;
	/** Probability that a polyploid variant is heterozygous. */
	private float hetRate=0.5f;

	/** Sum of inserted, deleted, substituted and junction bases for identity reporting. */
	private long mutationLengthAdded=0;
	/** Number of substitution events generated. */
	private long subsAdded=0;
	/** Number of insertion events generated. */
	private long insAdded=0;
	/** Number of deletion events generated. */
	private long delsAdded=0;
	/** Number of genome-fraction wraparound junctions generated. */
	private long junctionsAdded=0;

	/** If nonnegative, force mutation attempts after this many nonmutation loop iterations. */
	private int period=-1;

	/** Fraction of each input sequence retained before mutation. */
	private float genomeFraction=1;
	/** Sum of first-haplotype output lengths used for speed reporting. */
	private long basesRetained;

	/** Number of input reads processed. */
	private long readsProcessed=0;
	/** Number of input bases processed. */
	private long basesProcessed=0;

	/** Per-base substitution probability. */
	private float subRate=0;
	/** Per-base probability of attempting an indel event. */
	private float indelRate=0;
	/** Maximum sampled indel length. */
	private int maxIndel=1;
	/** Minimum reference spacing between random indel starts. */
	private int indelSpacing=10;
	/** If true, retry nucleotide indels that would create homopolymers. */
	private boolean banHomopolymers=false;
	/** Simple event probability used to decide whether to mutate the current base. */
	private final float errorRate;
	/** Approximate reported event-length-adjusted error rate. */
	private final float errorRate2;

	/** Seeded shared random source for variant generation. */
	private final Random randy;
	/** Seed passed to the shared thread-local random factory. */
	private long seed=-1;

	/*--------------------------------------------------------------*/

	public boolean errorState=false;
	/** Output overwrite flag. */
	private boolean overwrite=true;
	/** Output append flag. */
	private boolean append=false;

	/** Variant type constants matching {@link Var}. */
	private static final int SUB=Var.SUB, INS=Var.INS, DEL=Var.DEL;

	/** Diagnostic and status stream redirected by {@link PreParser}. */
	private java.io.PrintStream outstream=System.err;
	/** Verbose stream-loop logging flag. */
	public static boolean verbose=false;

}
