package tax;

import java.io.PrintStream;

import dna.Data;
import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.FastaReadInputStream;
import structures.ByteBuilder;
import tracker.ReadStats;

/**
 * Compresses and optimizes accession number storage and processing.
 * Processes NCBI accession-to-taxid mapping files to reduce size by removing
 * version numbers and optionally GI numbers while preserving taxonomic mappings.
 *
 * @author Brian Bushnell
 * @date April 4, 2017
 */
public class ShrinkAccession {
	
	public static void main(String[] args){
		Timer t=new Timer();
		ShrinkAccession x=new ShrinkAccession(args);
		x.process(t);
		
		//Close the print stream if it was redirected
		Shared.closeStream(x.outstream);
	}
	
	public ShrinkAccession(String[] args){
		
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, getClass(), false);
			args=pp.args;
			outstream=pp.outstream;
		}
		
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.setZipThreads(Shared.threads());
		if(Data.PIGZ()){
			ReadWrite.ZIPLEVEL=Tools.max(ReadWrite.ZIPLEVEL, 6);
		}
		
		Parser parser=new Parser();
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;

			if(parser.parse(arg, a, b)){
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
				ReadWrite.verbose=verbose;
			}else if(a.equals("seq") || a.equals("sequence")){
				mode=(Parse.parseBoolean(b) ? SEQUENCE : ASSEMBLY);
			}else if(a.equals("asm") || a.equals("assembly")){
				mode=(!Parse.parseBoolean(b) ? SEQUENCE : ASSEMBLY);
			}else if(a.equals("gi")){
				KEEP_GI_NUMBERS=Parse.parseBoolean(b);
			}else if(a.equals("outgi") || a.equals("giout") || a.equals("gi")){
				giOut=b;
			}else if(parser.in1==null && i==0 && Tools.looksLikeInputStream(arg)){
				parser.in1=arg;
			}else if(parser.out1==null && i==1 && !arg.contains("=")){
				parser.out1=arg;
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
				//				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		{//Process parser fields
			Parser.processQuality();
			
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;
			
			in=parser.in1;

			out=parser.out1;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in==null){throw new RuntimeException("Error - at least one input file is required.");}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2){
			ByteFile.FORCE_MODE_BF2=false;
			ByteFile.FORCE_MODE_BF1=true;
		}

		if(out!=null && out.equalsIgnoreCase("null")){out=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out)){
			outstream.println((out==null)+", "+out);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out+"\n");
		}

		ffout=FileFormat.testOutput(out, FileFormat.TXT, null, true, overwrite, append, false);
		ffoutGi=FileFormat.testOutput(giOut, FileFormat.TXT, null, true, overwrite, append, false);
		ffin=FileFormat.testInput(in, FileFormat.TXT, null, true, true);
		
	}
	
	void process(Timer t){
		
		ByteFile bf=ByteFile.makeByteFile(ffin);
		ByteStreamWriter bsw=new ByteStreamWriter(ffout);
		bsw.start();
		
		if(mode==SEQUENCE) {
			processSeq(bf, bsw);
		}else if(mode==ASSEMBLY){
			processAsm(bf, bsw);
		}else {
			throw new RuntimeException("Unknown mode: "+mode);
		}
		
		errorState|=bf.close();
		if(bsw!=null){errorState|=bsw.poisonAndWait();}
		
		t.stop();
		outstream.println("Discarded "+badLines+" lines.\n");
		outstream.println(Tools.timeLinesBytesProcessed(t, linesProcessed, charsProcessed, 8));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	void processSeq(ByteFile bf, ByteStreamWriter bsw){
		
		byte[] line=bf.nextLine();
		ByteBuilder bb=new ByteBuilder(10000);
		int columns=4;
		while(line!=null){
			if(Tools.startsWith(line, "accession\t")){
				bb.append(line);
				bb.nl();
			}else if(Tools.startsWith(line, "accession.version\ttaxid")){
				columns=2;
				bb.append("accession\t\ttaxid\t");//dummy header
				bb.nl();
			}else{
				charsProcessed+=line.length+1;
				linesProcessed++;
				
				final int tid=(columns==4 ? AccessionToTaxid.parseLineToTaxid(line, (byte)'\t') : 
					AccessionToTaxid.parseLineToTaxid_2col(line, (byte)'\t'));
				if(tid<1){
					badLines++;
				}else{
					int i=0;
					
					while(i<line.length){//Accession
						byte b=line[i];
						bb.append(b);
						i++;
						if(b=='\t'){break;}
					}
					
					if(columns==4){
						while(i<line.length){//Accession with decimal
							byte b=line[i];
							//						bb.append(b);
							i++;
							if(b=='\t'){break;}
						}
					}
					bb.append('\t');
					
					while(i<line.length){//Taxid
						byte b=line[i];
						bb.append(b);
						i++;
						if(b=='\t'){break;}
					}
					
					if(KEEP_GI_NUMBERS){
						if(line.length>i && Tools.isDigit(line[i])){//GI number or "na"
							while(i<line.length){
								byte b=line[i];
								bb.append(b);
								i++;
//								if(b=='\t'){break;}
							}
						}
					}
					bb.nl();
				}
			}
			if(bb.length()>8000){
				bsw.print(bb);
				bb.clear();
			}
			line=bf.nextLine();
		}
		if(bb.length()>0){
			bsw.print(bb);
			bb.clear();
		}
	}
	
//	col	name	sample
//	0	#assembly_accession	GCF_000001215.4
//	5	taxid	7227
//	23	assembly_type	haploid
//	25	genome_size	143706478
//	26	genome_size_ungapped	142553500
//	28	replicon_count	7
//	29	scaffold_count	1869
//	30	contig_count	2441
//	34	total_gene_count	17872
//	35	protein_coding_gene_count	13962
//	36	non_coding_gene_count	3543
	void processAsm(ByteFile bf, ByteStreamWriter bsw){
		
		byte[] line=bf.nextLine();
		ByteBuilder bb=new ByteBuilder(20000);
		LineParser1 lp=new LineParser1('\t');
//		int columns=4;
		while(line!=null && Tools.startsWith(line, "#")){
			if(Tools.startsWith(line, "#assembly_accession\t")){
				lp.set(line);
				bb.append(lp.parseByteArray(0)).tab();//accession
				bb.append(lp.parseByteArray(5)).tab();//tid
				bb.append(lp.parseByteArray(23)).tab();//ploidy?
				bb.append(lp.parseByteArray(25)).tab();//size
				bb.append(lp.parseByteArray(26)).tab();//size_ungapped
				bb.append(lp.parseByteArray(28)).tab();//chromosomes
				bb.append(lp.parseByteArray(29)).tab();//scafs
				bb.append(lp.parseByteArray(30)).tab();//contigs
				bb.append(lp.parseByteArray(34)).tab();//genes
				bb.append(lp.parseByteArray(35)).tab();//coding
				bb.append(lp.parseByteArray(36));//noncoding
				bb.nl();
			}
			charsProcessed+=line.length+1;
			linesProcessed++;
			line=bf.nextLine();
		}
		while(line!=null){
			charsProcessed+=line.length+1;
			linesProcessed++;

			if(Tools.startsWith(line, "#")){continue;}//Concatenated files or new header
			lp.set(line);

			bb.append(lp.parseByteArray(0));//accession
			bb.append(parseNum(lp, 5)).tab();//tid
			bb.append(lp.parseByteArray(23)).tab();//ploidy?
			bb.append(parseNum(lp, 25)).tab();//size
			bb.append(parseNum(lp, 26)).tab();//size_ungapped
			bb.append(parseNum(lp, 28)).tab();//chromosomes
			bb.append(parseNum(lp, 29)).tab();//scafs
			bb.append(parseNum(lp, 30)).tab();//contigs
			bb.append(parseNum(lp, 34)).tab();//genes
			bb.append(parseNum(lp, 35)).tab();//coding
			bb.append(parseNum(lp, 36));//noncoding
			bb.nl();
			if(bb.length()>16000){
				bsw.print(bb);
				bb.clear();
			}
			line=bf.nextLine();
		}
		if(bb.length()>0){
			bsw.print(bb);
			bb.clear();
		}
	}
	
	private long parseNum(LineParser1 lp, int field) {
		long x=-1;
		if(lp.termStartsWithLetter(field)) {return -1;}
//		try{
		x=lp.parseLong(field);
//		}catch(Throwable e){
//			e.printStackTrace();
//		}
		return x;
	}
	
	/*--------------------------------------------------------------*/
	
	
	/*--------------------------------------------------------------*/
	
	private String in=null;
	private String out=null;
	private String giOut=null;
	private int mode=SEQUENCE;

	long linesProcessed=0;
	long charsProcessed=0;
	long badLines=0;
	
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin;
	private final FileFormat ffout;
	private final FileFormat ffoutGi;
	
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public static boolean KEEP_GI_NUMBERS=true;
	public static final int SEQUENCE=1, ASSEMBLY=2;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
	
}
