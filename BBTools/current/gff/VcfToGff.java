package gff;

import java.io.PrintStream;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import shared.Parser;
import shared.PreParser;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import tracker.ReadStats;
import var2.VCFLine;

/**
 * Stripped out of GffLine into independent class.
 * @author Brian Bushnell
 * @date Sep 12, 2018
 *
 */
public class VcfToGff {

	/** Translates VCF to GFF */
	public static void main(String[] args){
		Timer t=new Timer();
		PrintStream outstream=System.err;
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, new Object() { }.getClass().getEnclosingClass(), false);
			args=pp.args;
			outstream=pp.outstream;
			t.outstream=outstream;
		}
		
		Parser parser=new Parser();
		String in=null;
		String out=null;
		boolean overwrite=true, append=false;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("in") || a.equals("vcf")){
				in=b;
			}else if(a.equals("out") || a.equals("gff")){
				out=b;
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(in==null && b==null && i==0 && Tools.canRead(arg)){
				in=arg;
			}else if(in==null && b==null && i==1){
				out=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		{//Process parser fields
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
		}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out+"\n");
		}

		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, in)){
			throw new RuntimeException("\nCan't read some input files.\n");  
		}

		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in, out)){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		translate(in, out, overwrite, append);
		t.stop("Time: \t");
	}
	
	/** Translates VCF to GFF */
	private static void translate(String in, String out, boolean overwrite, boolean append){
		//Create output FileFormat objects
		FileFormat ffout=FileFormat.testOutput(out, FileFormat.GFF, "gff", true, overwrite, append, false);

		//Create input FileFormat objects
		FileFormat ffin=FileFormat.testInput(in, FileFormat.VCF, "vcf", true, true);
		
		ByteFile bf=ByteFile.makeByteFile(ffin);
		ByteStreamWriter bsw=null;
		if(ffout!=null){
			bsw=new ByteStreamWriter(ffout);
			bsw.start();
		}
		
		ByteBuilder bb=new ByteBuilder(17000);
		bb.append("##gff-version 3\n");
		String header="#seqid	source	type	start	end	score	strand	phase	attributes";
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length>1){
				if(line[0]=='#'){
					if(Tools.startsWith(line, "##fileformat") || Tools.startsWith(line, "##FORMAT") || 
							Tools.startsWith(line, "##INFO") || Tools.startsWith(line, "#CHROM	POS")){
						//skip
					}else{
						int i=1;
						while(i<line.length && line[i]=='#'){i++;}
						i--;
						bb.append(line, i, line.length-i);
						bb.nl();
					}
				}else{
					if(header!=null){
						bb.append(header).append('\n');
						header=null;
					}
					VCFLine vline=new VCFLine(line);
					GffLine gline=new GffLine(vline);
					gline.appendTo(bb);
					bb.nl();
				}
			}
			if(bb.length()>=16384){
				if(bsw!=null){
					bsw.print(bb);
				}
				bb.clear();
			}
		}
		if(bb.length()>0){
			if(bsw!=null){
				bsw.print(bb);
			}
			bb.clear();
		}
		bf.close();
		if(bsw!=null){bsw.poisonAndWait();}
	}
	
}
