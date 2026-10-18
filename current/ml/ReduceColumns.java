package ml;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import parse.LineParser1;
import shared.KillSwitch;
import shared.Timer;
import shared.Tools;
import structures.IntList;
import structures.ListNum;

/**
 * Rewrites tab-delimited matrix rows with only selected columns.
 * <p>
 * Header lines beginning with {@code '#'} are skipped.  The first data row
 * resolves the last supplied open-ended {@code N+} column range and appends
 * that expansion after all explicit selections, then emits a BBTools
 * {@code #dims} header where the last selected column is the single output
 * column and every earlier selected column is treated as input.
 *
 * @author Brian Bushnell, Chloe
 * @date June 3, 2025
 */
public class ReduceColumns{

	/**
	 * Parses {@code in out colspec...}, then streams selected columns to the
	 * destination.  Column specs may be individual zero-based indexes or closed
	 * inclusive ranges such as {@code 0-8}; they preserve the listed order and
	 * duplicates.  At most one open range is retained: the last {@code N+} spec
	 * expands against the first data row's column count after all explicit
	 * selections.  Every selected index must exist in every data row.
	 *
	 * @param args input path, output path, and one or more column specs
	 */
	public static void main(String[] args){

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
		final boolean writeError=bsw.poisonAndWait();
		assert(!writeError) : KillSwitch.assertDie("Could not write reduced columns output: "+out);

		t.stop();
		System.err.println(Tools.timeLinesBytesProcessed(t, linesIn, bytesIn, 8));
	}

}
