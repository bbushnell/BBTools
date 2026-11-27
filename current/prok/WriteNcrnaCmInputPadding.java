package prok;

import java.io.File;
import java.io.PrintStream;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Timer;

/** Writes explicitly fixed input-padding resources from one consensus read.
 * Uses the basic template's argument/lifecycle shape and the existing native
 * sequence loader. No teaching calls, fitted support or payload hashes enter.
 * @author Brian Bushnell, Raiden
 */
public final class WriteNcrnaCmInputPadding{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, WriteNcrnaCmInputPadding.class, false);
		try{
			new WriteNcrnaCmInputPadding(pp.args, pp.outstream).process(timer);
		}finally{Shared.closeStream(pp.outstream);}
	}

	private WriteNcrnaCmInputPadding(final String[] args, final PrintStream outstream_){
		outstream=outstream_;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("pad5")){
				pad5=pads(b);
			}else if(a.equals("pad3")){
				pad3=pads(b);
			}else if(!parser.parse(arg, a, b)){
				throw new IllegalArgumentException("Unknown fixed-padding export option: "+arg);
			}
		}
		require(parser.in1!=null && parser.in2==null && parser.out1!=null && parser.out2==null
			&& !parser.append && !parser.overwrite && parser.maxReads<0,
			"Fixed-padding export requires one complete consensus FASTA and a fresh output directory; no append, overwrite, paired input or read limit");
		require(pad5!=null && pad3!=null && pad5.length==pad3.length && pad5.length>0,
			"pad5= and pad3= are parallel nonempty lists of fixed sense-oriented base counts");
		for(int i=0; i<pad5.length; i++){
			for(int j=0; j<i; j++){
				require(pad5[i]!=pad5[j] || pad3[i]!=pad3[j], "Repeated fixed-padding pair: "+pad5[i]+","+pad3[i]);
			}
		}
		input=parser.in1;
		output=new File(parser.out1);
		require(!output.exists(), "Preserve existing fixed-padding exports: "+output);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Export             ----------------*/
	/*--------------------------------------------------------------*/

	private void process(final Timer timer){
		assert(input!=null && output!=null && pad5.length==pad3.length) :
			"Validated parallel pad pairs and one consensus axis define every output resource";
		final NcrnaSupport.Library library=NcrnaSupport.loadLibrary(input);
		require(output.mkdir(), "Could not create the fresh fixed-padding directory: "+output);
		for(int i=0; i<pad5.length; i++){
			final File file=new File(output, "padding_"+pad5[i]+"_"+pad3[i]+".tsv");
			NcrnaCmInputPadding.writeFixed(file.getPath(), pad5[i], pad3[i], library.names, library.bases);
		}
		timer.stop();
		outstream.println("Models: "+library.names.length+"; fixed padding pairs: "+pad5.length+"; time: "+timer);
		outstream.println("NCRNA_CM_FIXED_PADDING_EXPORT_COMPLETE");
	}

	private static int[] pads(final String value){
		require(value!=null && !value.isEmpty(), "Fixed padding needs a nonempty comma-separated list");
		final String[] terms=value.split(",", -1);
		final int[] result=new int[terms.length];
		for(int i=0; i<result.length; i++){
			require(!terms[i].isEmpty(), "Empty fixed-padding list entry");
			final long bases=Parse.parseKMG(terms[i]);
			require(bases>=0 && bases<=Integer.MAX_VALUE, "Fixed padding must fit a nonnegative int base count: "+terms[i]);
			result[i]=(int)bases;
		}
		return result;
	}

	private static void require(final boolean condition, final String why){
		if(!condition){throw new IllegalArgumentException(why);}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final String input;
	private final File output;
	private final PrintStream outstream;
	private int[] pad5, pad3;
}
