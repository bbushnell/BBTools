package var2;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import fileIO.ByteFile;
import fileIO.ByteFile1;
import fileIO.ByteFile2;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Inverts a VCF file produced by MutateGenome (mutate.sh).
 * Swaps ref/alt alleles, flips INS/DEL types, and adjusts coordinates
 * from original-genome space to mutant-genome space.
 *
 * This allows comparison of variant calls made by mapping real reads
 * to a mutant genome against the known ground truth of the mutations.
 *
 * @author Brian Bushnell, Chloe
 * @date June 2026
 */
public class InvertVCF {

	public static void main(String[] args){
		Timer t=new Timer();
		InvertVCF x=new InvertVCF(args);
		x.process(t);
		Shared.closeStream(x.outstream);
	}

	public InvertVCF(String[] args){

		{
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

			if(a.equals("verbose")){
				verbose=parse.Parse.parseBoolean(b);
				ByteFile1.verbose=verbose;
				ByteFile2.verbose=verbose;
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}

		{
			in1=parser.in1;
			out1=parser.out1;
			overwrite=parser.overwrite;
			append=parser.append;
		}

		if(in1==null){throw new RuntimeException("Error - an input file is required.");}

		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}

		if(!Tools.testInputFiles(false, true, in1)){
			throw new RuntimeException("\nCan't read input file.\n");
		}

		ffin1=FileFormat.testInput(in1, FileFormat.TXT, null, true, true);
		ffout1=FileFormat.testOutput(out1, FileFormat.TXT, null, true, overwrite, append, false);
	}

	void process(Timer t){

		ByteFile bf=ByteFile.makeByteFile(ffin1);

		ArrayList<byte[]> headerLines=new ArrayList<byte[]>();
		ArrayList<byte[]> dataLines=new ArrayList<byte[]>();

		byte[] line=bf.nextLine();
		while(line!=null){
			if(line.length>0){
				linesIn++;
				bytesIn+=line.length;
				if(line[0]=='#'){
					headerLines.add(line);
				}else{
					dataLines.add(line);
				}
			}
			line=bf.nextLine();
		}
		errorState|=bf.close();

		LinkedHashMap<String, Long> netShiftPerScaf=new LinkedHashMap<String, Long>();
		for(byte[] dline : dataLines){
			String[] fields=new String(dline).split("\t");
			String chrom=fields[0];
			String refStr=fields[3];
			String altStr=fields[4];
			int refLen=refStr.length();
			int altLen=altStr.length();
			int delta=altLen-refLen;
			Long prev=netShiftPerScaf.get(chrom);
			netShiftPerScaf.put(chrom, (prev==null ? 0 : prev)+delta);
		}

		ByteStreamWriter bsw=null;
		if(ffout1!=null){
			bsw=new ByteStreamWriter(ffout1);
			bsw.start();
		}

		ByteBuilder bb=new ByteBuilder();
		for(byte[] hline : headerLines){
			String hs=new String(hline);
			if(hs.startsWith("##contig=<ID=")){
				String adjusted=adjustContigLength(hs, netShiftPerScaf);
				bb.append(adjusted).append('\n');
			}else if(hs.startsWith("##Program=")){
				bb.append(hs).append('\n');
				bb.append("##InvertedBy=InvertVCF").append('\n');
			}else{
				bb.append(hline).append('\n');
			}
			headerLinesOut++;
		}
		if(bsw!=null && bb.length()>0){
			bsw.print(bb);
			bb.clear();
		}

		String prevChrom=null;
		long cumShift=0;

		for(byte[] dline : dataLines){
			String[] fields=new String(dline).split("\t");
			String chrom=fields[0];
			int pos=Integer.parseInt(fields[1]);
			String id_f=fields[2];
			String refStr=fields[3];
			String altStr=fields[4];

			if(!chrom.equals(prevChrom)){
				cumShift=0;
				prevChrom=chrom;
			}

			int newPos=pos+(int)cumShift;

			int refLen=refStr.length();
			int altLen=altStr.length();
			int delta=altLen-refLen;

			String newRef=altStr;
			String newAlt=refStr;

			String infoStr=fields[7];
			String newInfo=invertInfo(infoStr, cumShift);

			cumShift+=delta;

			bb.append(chrom).append('\t');
			bb.append(newPos).append('\t');
			bb.append(id_f).append('\t');
			bb.append(newRef).append('\t');
			bb.append(newAlt).append('\t');
			for(int i=5; i<fields.length; i++){
				if(i==7){
					bb.append(newInfo);
				}else{
					bb.append(fields[i]);
				}
				if(i<fields.length-1){bb.append('\t');}
			}
			bb.append('\n');
			variantLinesOut++;

			if(bb.length()>=64000){
				if(bsw!=null){bsw.print(bb);}
				bb.clear();
			}
		}

		if(bb.length()>0 && bsw!=null){
			bsw.print(bb);
			bb.clear();
		}

		if(bsw!=null){errorState|=bsw.poisonAndWait();}

		t.stop();
		outstream.println(Tools.timeLinesBytesProcessed(t, linesIn, bytesIn, 8));
		outstream.println();
		outstream.println("Header Lines Out:  \t"+headerLinesOut);
		outstream.println("Variant Lines Out: \t"+variantLinesOut);

		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}

	/**
	 * Adjusts the length in a ##contig header line by the net indel shift
	 * for that scaffold.
	 */
	private String adjustContigLength(String headerLine, LinkedHashMap<String, Long> netShiftPerScaf){
		int idStart=headerLine.indexOf("ID=")+3;
		int idEnd=headerLine.indexOf(',', idStart);
		if(idEnd<0){idEnd=headerLine.indexOf('>', idStart);}
		String scafName=headerLine.substring(idStart, idEnd);

		int lenStart=headerLine.indexOf("length=");
		if(lenStart<0){return headerLine;}
		lenStart+=7;
		int lenEnd=lenStart;
		while(lenEnd<headerLine.length() && Character.isDigit(headerLine.charAt(lenEnd))){lenEnd++;}
		long oldLen=Long.parseLong(headerLine.substring(lenStart, lenEnd));

		Long shift=netShiftPerScaf.get(scafName);
		if(shift==null){shift=0L;}
		long newLen=oldLen+shift;

		return headerLine.substring(0, lenStart)+newLen+headerLine.substring(lenEnd);
	}

	/**
	 * Inverts the INFO field: updates STA, STO coordinates with cumShift,
	 * and flips TYP between INS and DEL.
	 */
	private String invertInfo(String info, long cumShift){
		String[] parts=info.split(";");
		StringBuilder sb=new StringBuilder();
		for(int i=0; i<parts.length; i++){
			if(i>0){sb.append(';');}
			String part=parts[i];
			if(part.startsWith("STA=")){
				long oldSta=Long.parseLong(part.substring(4));
				sb.append("STA=").append(oldSta+cumShift);
			}else if(part.startsWith("STO=")){
				long oldSto=Long.parseLong(part.substring(4));
				sb.append("STO=").append(oldSto+cumShift);
			}else if(part.startsWith("TYP=")){
				String oldType=part.substring(4);
				if(oldType.equals("INS")){
					sb.append("TYP=DEL");
				}else if(oldType.equals("DEL")){
					sb.append("TYP=INS");
				}else{
					sb.append(part);
				}
			}else{
				sb.append(part);
			}
		}
		return sb.toString();
	}

	private String in1=null;
	private String out1=null;

	private final FileFormat ffin1;
	private final FileFormat ffout1;

	private long linesIn=0;
	private long bytesIn=0;
	private long headerLinesOut=0;
	private long variantLinesOut=0;

	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=true;
	private boolean append=false;
}
