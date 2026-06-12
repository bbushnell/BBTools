package ml;

import java.io.PrintStream;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import parse.PreParser;
import shared.Tools;
import structures.ByteBuilder;

/**
 * Reads a neural network (.bbnet) and writes it back out, purely to convert its
 * on-disk format.  The output coding follows {@link CellNet#codingA48Out} (A48 by
 * default), so an old decimal network becomes A48 -- lossless for the stored
 * weights and roughly 30% smaller.  The weights themselves are never altered;
 * this is a format round-trip, not a re-training.
 *
 * Usage:  netconvert.sh in=&lt;old.bbnet&gt; out=&lt;new.bbnet&gt;
 *
 * @author UMP45
 */
public class NetConvert {

	public static void main(String[] args){
		PreParser pp=new PreParser(args, NetConvert.class, false);
		args=pp.args;
		PrintStream outstream=pp.outstream;

		String in=null, out=null;
		boolean overwrite=true;
		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(a.equals("in") || a.equals("net") || a.equals("netin")){in=b;}
			else if(a.equals("out") || a.equals("netout")){out=b;}
			else if(a.equals("a48")){CellNet.codingA48Out=Parse.parseBoolean(b);}
			else if(a.equals("overwrite") || a.equals("ow")){overwrite=Parse.parseBoolean(b);}
			else{outstream.println("Unknown parameter "+arg); assert(false) : "Unknown parameter "+arg;}
		}
		if(in==null || out==null){
			throw new RuntimeException("Usage: netconvert.sh in=<old.bbnet> out=<new.bbnet>");
		}
		if(!Tools.testOutputFiles(overwrite, false, false, out)){
			throw new RuntimeException("Can't write to "+out+" (overwrite="+overwrite+").");
		}

		// Read (coding auto-detected from the #coding header, or decimal for old headerless nets).
		CellNet net=CellNetParser.load(in);
		assert(net!=null) : "Failed to load network: "+in;

		// Write (coding = CellNet.codingA48Out, A48 by default).
		FileFormat ffout=FileFormat.testOutput(out, FileFormat.BBNET, null, true, overwrite, false, false);
		ByteStreamWriter bsw=new ByteStreamWriter(ffout);
		bsw.start();
		ByteBuilder bb=net.toBytes();
		bsw.println(bb);
		bsw.poisonAndWait();

		outstream.println("Converted "+in+" -> "+out
				+"  (coding="+(CellNet.codingA48Out ? "A48" : "decimal")+", "+bb.length()+" bytes)");
	}
}
