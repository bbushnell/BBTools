package ml;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import shared.LineParser1;
import shared.Timer;
import shared.Tools;
import structures.IntList;
import structures.ListNum;

public class ReduceColumns {

	public static void main(String[] args) {

		Timer t=new Timer();
		String in=args[0];
		String out=args[1];
		IntList columns=new IntList();
		for(int i=2; i<args.length; i++) {
			columns.add(Integer.parseInt(args[i]));
		}
		columns.shrink();
		
		ByteFile bf=ByteFile.makeByteFile(in, true);
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(out, true, false, true);
		
		LineParser1 lp=new LineParser1('\t');
		boolean header=false;
		
		bsw.print("#dims").tab().print(columns.size()-1).tab().print(1).nl();
		header=true;

		long linesIn=0, bytesIn=0;
		for(ListNum<byte[]> ln=bf.nextList(); ln!=null; ln=bf.nextList()) {
			for(byte[] line : ln) {
				linesIn++;
				bytesIn+=line.length;
				lp.set(line);
				if(lp.startsWith('#')) {
//					int a=lp.parseInt(1);
//					int b=lp.parseInt(2);//TODO: figure out new dims
//					if(!header && lp.startsWith("#dims\t")) {
//						bsw.println("#dims").tab().print(columns.size()-1).tab().print(1);
//						header=true;
//					}
				}else {
					for(int i=0; i<columns.size; i++) {
						if(i>0) {bsw.tab();}
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
