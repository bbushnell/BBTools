package gff;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import shared.LineParser1;
import shared.Tools;
import structures.ByteBuilder;

public class GtfLine {
	
	public static void main(String[] args) {
		String in=args[0];
		String out=args.length>1 ? args[1] : "stdout.gff";
		ByteFile bf=ByteFile.makeByteFile(in, true);
		ByteStreamWriter bsw=ByteStreamWriter.makeBSW(out, true, false, true);
		
		ByteBuilder bb=new ByteBuilder();
		bb.append("##gff-version 3\n");
		bb.append("#seqid	source	type	start	end	score	strand	phase	attributes\n");
		bsw.print(bb);
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()) {
			if(Tools.startsWith(line, '#')) {continue;}
			GtfLine gtf=new GtfLine(line);
			GffLine gff=new GffLine(gtf);
			gff.appendTo(bb.clear());
			bsw.print(bb.nl());
		}
		bsw.poisonAndWait();
	}
	
	public GtfLine(byte[] line){
		
		LineParser1 lp=new LineParser1('\t');
		lp.set(line);
		
	    seqname=lp.parseString(0);
	    source=lp.parseString(1);
	    feature=lp.parseString(2);
	    start=lp.parseInt(3);
	    end=lp.parseInt(4);
	    score=(lp.termEquals('.', 5) ? -1 : lp.parseFloat(5));
	    strand=lp.parseByte(6, 0);
	    frame=(lp.termEquals('.', 7) ? -1 : lp.parseInt(7));
	    attribute=lp.parseString(8);
//		if(start==267921) {System.err.println(lp+" -> "+attribute);}
	}
	
    String seqname;
    String source;
    String feature;
    int start;
    int end;
    float score;
    byte strand;
    int frame;
    String attribute;
	
}
