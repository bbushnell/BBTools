package ddl;

import java.io.PrintStream;
import java.util.ArrayList;

import shared.Timer;

/**
 * Converts legacy DDL sketch files (absolute encoding, no #offset header)
 * to the new format (relative encoding with #offset header).
 *
 * The conversion is lossless: DDLLoader.loadFile reads absolute values via
 * fromArray(offset=-1) which auto-converts to relative encoding internally.
 * DDLLoader.writeFile then writes #offset + relative values.
 *
 * Usage: DDLConverter in=refseq_ddls.tsv.gz out=refseq_ddls_v2.tsv.gz [k=31]
 *
 * @author Brian Bushnell, UMP45
 * @date April 23, 2026
 */
public class DDLConverter {

	public static void main(String[] args){
		String in=null, out=null;
		int k=31;
		boolean overwrite=false;

		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(a.equals("in")){in=b;}
			else if(a.equals("out")){out=b;}
			else if(a.equals("k")){k=Integer.parseInt(b);}
			else if(a.equals("overwrite") || a.equals("ow")){overwrite=true;}
			else{
				outstream.println("Unknown parameter: "+arg);
				assert(false) : "Unknown parameter: "+arg;
			}
		}

		if(in==null || out==null){
			outstream.println("Usage: DDLConverter in=<input.tsv[.gz]> out=<output.tsv[.gz]> [k=31] [overwrite]");
			System.exit(1);
		}

		Timer t=new Timer();
		outstream.println("Loading "+in+" ...");
		ArrayList<DDLRecord> records=DDLLoader.loadFile(in, k);
		outstream.println("Loaded "+records.size()+" records.");

		outstream.println("Writing "+out+" ...");
		DDLLoader.writeFile(records, out, overwrite, k, -1);
		t.stop();

		outstream.println("Converted "+records.size()+" records to relative encoding with #offset header.");
		outstream.println("Time: \t"+t);
	}

	private static final PrintStream outstream=System.err;
}
