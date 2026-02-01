package gff;

import java.io.PrintStream;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parser;
import parse.PreParser;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;
import tracker.ReadStats;
import var2.VCFLine;

/**
 * Converts Variant Call Format (VCF) files to Gene Feature Format (GFF) files.
 * A command-line utility that parses input VCF files and transforms variant
 * information into GFF format, preserving metadata and generating a standard
 * GFF3 output with genomic feature annotations.
 *
 * Key features:
 * - Translates VCF variant records to GFF3 format
 * - Preserves VCF metadata by selectively copying header comments
 * - Supports input/output file specification via command-line arguments
 * - Generates a standard GFF3 header with version 3
 * - Implements robust file handling with overwrite and append options
 * - Uses byte-level file processing for memory efficiency
 *
 * @author Brian Bushnell
 */
public class VcfToGff {

	/**
	 * Main entry point for VCF to GFF conversion.
	 * Parses command-line arguments, validates input/output files, and initiates
	 * the translation process. Supports flexible argument parsing with input
	 * and output file specification.
	 *
	 * @param args Command-line arguments including:
	 * in/vcf=input.vcf - Input VCF file
	 * out/gff=output.gff - Output GFF file
	 */
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
	
	/**
	 * Performs the actual VCF to GFF conversion process.
	 * Reads VCF file line-by-line, converts VCFLine objects to GffLine objects,
	 * and writes the result to the output file. Preserves relevant VCF header
	 * information while filtering out format-specific headers.
	 *
	 * @param in Input VCF file path
	 * @param out Output GFF file path
	 * @param overwrite Whether to overwrite existing output files
	 * @param append Whether to append to existing output files
	 */
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
