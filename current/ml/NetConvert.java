package ml;

import java.io.PrintStream;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import parse.PreParser;
import shared.KillSwitch;
import shared.Tools;
import structures.ByteBuilder;

/** Converts a .bbnet file to the currently selected network writer format.
 * Input coding is autodetected from {@code #coding} headers, or falls back to
 * decimal parsing for legacy headerless networks.  Output uses the process-wide
 * {@link CellNet} writer flags, defaulting to lossless A48 unless they were
 * changed earlier in the same JVM or overridden by this command.  Writer-flag
 * parsing mutates those globals; the weights themselves are not retrained or
 * otherwise transformed.
 * <p>Usage: {@code netconvert.sh in=<old.bbnet> out=<new.bbnet>}
 *
 * @author UMP45
 */
public class NetConvert{

	/** Parses writer-format flags, loads one network and writes it in the selected format. */
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
			else if(a.equals("outsparse")){CellNet.OUT_SPARSE=Parse.parseBoolean(b); if(CellNet.OUT_SPARSE){CellNet.OUT_DENSE=false;}}
			else if(a.equals("outdense")){CellNet.OUT_DENSE=Parse.parseBoolean(b); if(CellNet.OUT_DENSE){CellNet.OUT_SPARSE=false;}}
			else if(a.equals("outhex")){CellNet.OUT_HEX=Parse.parseBoolean(b); if(CellNet.OUT_HEX){CellNet.OUT_DELTA_A48=false;}}
			else if(a.equals("outdelta") || a.equals("outdeltaa48") || a.equals("deltaa48")){
				CellNet.OUT_DELTA_A48=Parse.parseBoolean(b);
				if(CellNet.OUT_DELTA_A48){CellNet.OUT_HEX=false;}
			}
			else if(a.equals("overwrite") || a.equals("ow")){overwrite=Parse.parseBoolean(b);}
			else{outstream.println("Unknown parameter "+arg); assert(false) : "Unknown parameter "+arg;}
		}
		if(in==null || out==null){
			throw new RuntimeException("Usage: netconvert.sh in=<old.bbnet> out=<new.bbnet>");
		}
		if(!Tools.testOutputFiles(overwrite, false, false, out)){
			throw new RuntimeException("Can't write to "+out+" (overwrite="+overwrite+").");
		}

		// Input coding auto-detected from #coding, or decimal for old headerless nets.
		CellNet net=CellNetParser.load(in);
		assert(net!=null) : "Failed to load network: "+in;

		// Output coding and dense/sparse/index flags are CellNet's current static writer policy.
		FileFormat ffout=FileFormat.testOutput(out, FileFormat.BBNET, null, true, overwrite, false, false);
		ByteStreamWriter bsw=new ByteStreamWriter(ffout);
		bsw.start();
		ByteBuilder bb=net.toBytes();
		bsw.println(bb);
		final boolean writeError=bsw.poisonAndWait();
		assert(!writeError) : KillSwitch.assertDie("Could not write converted network output: "+out);

		outstream.println("Converted "+in+" -> "+out
				+"  (coding="+(CellNet.codingA48Out ? "A48" : "decimal")+", "+bb.length()+" bytes)");
	}
}
