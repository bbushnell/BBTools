package synth;

import java.io.PrintStream;
import java.util.Random;

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
import structures.ByteBuilder;

/**
 * @author Brian Bushnell
 * @date Jan 3, 2013
 *
 */
public class RandomGenome {
	
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		RandomGenome x=new RandomGenome(args);
		
		//Run the object
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public RandomGenome(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b!=null && b.equalsIgnoreCase("null")){b=null;}

			if(a.equals("chroms")){
				chroms=Parse.parseIntKMG(b);
			}else if(a.equals("len") || a.equals("length") || a.equals("size")){
				totalLength=Parse.parseKMG(b);
			}else if(a.equals("pad")){
				pad=Tools.max(0, Parse.parseIntKMG(b));
			}else if(a.equals("gc")){
				gc=Float.parseFloat(b);
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ReadWrite.verbose=verbose;
			}else if(a.equals("nohomopolymers") || a.equals("banhomopolymers") || a.equals("nopoly")){
				noPoly=Parse.parseBoolean(b);
			}else if(a.equals("includestop") || a.equals("stop")){
				includeStop=Parse.parseBoolean(b);
			}else if(a.equals("seed")){
				seed=Long.parseLong(b);
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			overwrite=parser.overwrite;
			append=parser.append;

			out=parser.out1;
		}

		wrap=Shared.FASTA_WRAP;
		assert(wrap>0) : "Wrap is too small.";
		assert(chroms>0) : "Chroms must be greater than 0.";
		assert(totalLength>=chroms) : "Length must be at least chroms.";
		assert(2*pad+totalLength/chroms<Shared.MAX_ARRAY_LEN) : "Length per chrom must be at most "+Shared.MAX_ARRAY_LEN;
		chromLength=(int)(totalLength/chroms);

		if(out!=null && out.equalsIgnoreCase("null")){out=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out+"\n");
		}

		ffout=FileFormat.testOutput(out, FileFormat.FA, null, true, overwrite, append, false);
		
		randy=Shared.threadLocalRandom(seed);
	}
	
	void process(Timer t){
		if(Shared.AMINO_IN){
			processAmino(t);
		}else{
			processNucleotide(t);
		}
	}
	
	void processNucleotide(Timer t){
		
		ByteStreamWriter bsw=new ByteStreamWriter(ffout);
		bsw.start();
		
		for(int chrom=1; chrom<=chroms; chrom++){
			bsw.print('>').print("contig").println(chrom);
			ByteBuilder bb=new ByteBuilder(wrap+1);
			byte prev='N';
			final int max=chromLength+2*pad;
			final int pad2=chromLength+pad;
			if(gc==0.5f){
				for(int i=0; i<max; ){
					for(int j=0; j<wrap && i<max; i++, j++){
						byte b;
						if(i<pad || i>=pad2){b='N';}
						else{
							b=AminoAcid.numberToBase[randy.nextInt(4)];
							while(noPoly && b==prev){b=AminoAcid.numberToBase[randy.nextInt(4)];}
						}
						bb.append(b);
						prev=b;
					}
					bb.nl();
					bsw.print(bb);
					bb.clear();
				}
			}else{
				for(int i=0; i<max; ){
					for(int j=0; j<wrap && i<max; i++, j++){
						boolean at=randy.nextFloat()>=gc;
						char b;
						if(i<pad || i>=pad2){b='N';}
						else{
							if(at){
								b=randy.nextBoolean() ? 'A' : 'T';
							}else{
								b=randy.nextBoolean() ? 'C' : 'G';
							}
							while(noPoly && b==prev){
								if(at){
									b=randy.nextBoolean() ? 'A' : 'T';
								}else{
									b=randy.nextBoolean() ? 'C' : 'G';
								}
							}
						}
						bb.append(b);
						prev=(byte)b;
					}
					bb.nl();
					bsw.print(bb);
					bb.clear();
				}
			}
		}
		bsw.poisonAndWait();
	}
	
	void processAmino(Timer t){
		
		ByteStreamWriter bsw=new ByteStreamWriter(ffout);
		bsw.start();

		final byte[] acids=AminoAcid.numberToAcid;
		final int limit=(includeStop ? acids.length : acids.length-1);
		for(int chrom=1; chrom<=chroms; chrom++){
			bsw.print('>').print("gene").println(chrom);
			ByteBuilder bb=new ByteBuilder(wrap+1);
			byte prev='X';
			final int max=chromLength+2*pad;
			final int pad2=chromLength+pad;
			for(int i=0; i<max; ){
				for(int j=0; j<wrap && i<max; i++, j++){
					byte b;
					if(i<pad || i>=pad2){b='X';}
					else{
						b=acids[randy.nextInt(limit)];
						while(noPoly && b==prev){b=acids[randy.nextInt(limit)];}
					}
					bb.append(b);
					prev=b;
				}
				bb.nl();
				bsw.print(bb);
				bb.clear();
			}
		}
		bsw.poisonAndWait();
	}
	
	/*--------------------------------------------------------------*/
	
	private String out=null;
	
	int chroms=1;
	long totalLength=1000000;
	float gc=0.5f;
	final int chromLength;
	final int wrap;
	int pad=0;
	boolean noPoly=false;
	boolean includeStop=false;
	long seed=-1;
	
	/*--------------------------------------------------------------*/

	final Random randy;
	
	private long linesOut=0;
	private long bytesOut=0;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffout;
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
