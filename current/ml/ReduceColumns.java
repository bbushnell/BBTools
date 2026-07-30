package ml;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import shared.Timer;
import shared.Tools;
import structures.IntList;
import structures.ListNum;

/**
 * Command-line utility for extracting specific columns from tab-delimited text files
 * with efficient memory management. Processes input files by selectively extracting
 * user-specified columns, writing results to a new output file with low memory overhead.
 * Useful for column extraction from large tabular data files and data preprocessing
 * for machine learning datasets.
 *
 * @author Brian Bushnell, Chloe
 * @date June 3, 2025
 */
public class ReduceColumns {

	/**
	 * Main entry point for the column reduction utility.
	 * Expects input file path, output file path, followed by column specs to extract.
	 * Supports individual columns (5), ranges (0-8), and open-ended ranges (17+).
	 *
	 * @param args Command line arguments: [input_file] [output_file] [colspec] [colspec] ...
	 */
	public static void main(String[] args) {

		Timer t=new Timer();
		String in=args[0];
		String out=args[1];
		IntList columns=new IntList();
		int openRangeStart=-1;
		for(int i=2; i<args.length; i++){
			String arg=args[i];
			if(arg.endsWith("+")){
				openRangeStart=Integer.parseInt(arg.substring(0, arg.length()-1));
			}else if(arg.contains("-")){
				int hyphen=arg.indexOf('-');
				int from=Integer.parseInt(arg.substring(0, hyphen));
				int to=Integer.parseInt(arg.substring(hyphen+1));
				for(int c=from; c<=to; c++){columns.add(c);}
			}else{
				columns.add(Integer.parseInt(arg));
			}
		}

		ByteFile bf=ByteFile.makeByteFile(in, true);
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(out, true, false, true);

		LineParser1 lp=new LineParser1('\t');
		boolean header=false;

		long linesIn=0, bytesIn=0;
		for(ListNum<byte[]> ln=bf.nextList(); ln!=null; ln=bf.nextList()){
			for(byte[] line : ln){
				linesIn++;
				bytesIn+=line.length;
				lp.set(line);
				if(lp.startsWith('#')){
					//Skip input headers; we write our own after resolving open ranges
				}else{
					if(!header){
						if(openRangeStart>=0){
							int totalCols=lp.terms();
							for(int c=openRangeStart; c<totalCols; c++){columns.add(c);}
							openRangeStart=-1;
						}
						columns.shrink();
						bsw.print("#dims").tab().print(columns.size()-1).tab().print(1).nl();
						header=true;
					}
					for(int i=0; i<columns.size; i++){
						if(i>0){bsw.tab();}
						bsw.print(lp.parseByteArray(columns.get(i)));
					}
					bsw.println();
				}
			}
		}
		bsw.poisonAndWait();

		t.stop();
		System.err.println(Tools.timeLinesBytesProcessed(t, linesIn, bytesIn, 8));
	}

}
