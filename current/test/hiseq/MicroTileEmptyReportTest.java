package test.hiseq;

import java.util.Locale;

import hiseq.MicroTile;
import structures.ByteBuilder;

/** Literal empty, SAM-only and nonempty reporting controls for HIL-109.
 * @author Brian Bushnell, Sucrose */
public class MicroTileEmptyReportTest{

	/** Checks baseline/fixed zero-base inference with/without a fit, plus a nonempty control and formatter probe. */
	public static void main(String[] args){
		boolean fixed=false;
		for(String arg : args){
			if(arg.equals("phase=fixed")){fixed=true;}
			else if(!arg.equals("phase=baseline")){throw new IllegalArgumentException(arg);}
		}
		Locale.setDefault(Locale.ROOT);
		MicroTile.TRACK_CYCLES=false;
		MicroTile.shortHeader=false;
		//Diagnostic reference expression, not a capture of MicroTile's local E2 variable.
		final double emptyReference=1-Math.pow(1.0, Double.POSITIVE_INFINITY);
		final String formattedNaN=new ByteBuilder().append(emptyReference, 5).toString();
		System.out.println("reference_empty_arithmetic="+emptyReference+"\tfixed_decimal_text="+formattedNaN);
		require(Double.isNaN(emptyReference), "The zero-base reference expression must expose pow(1,Infinity)");
		require("0.00000".equals(formattedNaN), "ByteBuilder's fixed-decimal path casts NaN to zero digits: "+formattedNaN);
		final String[] names=MicroTile.header().split("\t");
		require("inferredErrorRate2".equals(names[24]), "Report column24 must be inferredErrorRate2");
		System.out.println("case\treads\tbases\taligned_reads\taligned_bases\tinferred_error2\tresult");
		for(int scenario=0; scenario<4; scenario++){
			final MicroTile tile=new MicroTile(1, 1101, 0, 499, 0, 499);
			if(scenario==1){
				//AnalyzeFlowCell.processSamLine increments alignment counts without read/base totals.
				tile.alignedReadCount=1; tile.alignedBaseCount=4;
			}else if(scenario==2){
				tile.readCount=1; tile.baseCount=4; tile.acgtn[0]=4;
			}
			final String name=scenario==0 ? "empty" : scenario==1 ? "sam-only" :
				scenario==2 ? "nonempty-control" : "empty-with-fit";
			//For zero bases and a positive fitted rate, pow(0.99,Infinity)=0, so the old E2 is1.
			final double[] fit=scenario==3 ? new double[] {0.01, 0} : null;
			final String text=tile.toText(new ByteBuilder(), 31, 0, fit, null).toString().trim();
			final String[] row=text.split("\t");
			System.out.println("OBSERVED\t"+name+"\trow="+text);
			System.out.println("OBSERVED\t"+name+"\tinferredErrorRate2="+(row.length>24 ? row[24] : "MISSING"));
			require(row.length>=25, "Core report fields must be present");
			require(row[6].equals(scenario==2 ? "1" : "0") && row[7].equals(scenario==2 ? "4" : "0"),
				"Literal input read/base counts changed");
			require(row[8].equals(scenario==1 ? "1" : "0") && row[9].equals(scenario==1 ? "4" : "0"),
				"Literal alignment counts changed");
			final double inferred=Double.parseDouble(row[24]);
			final double expected=!fixed && scenario==3 ? 1 : 0;
			require(inferred==expected, "Literal zero-base/nonempty report oracle: "+name+"="+row[24]+", expected="+expected);
			if(!fixed){require(row[24].equals(scenario<2 ? "0.00000" : scenario==3 ? "1" : "0"),
				"Baseline must distinguish formatted NaN digits from the integer-zero fast path: "+name+"="+row[24]);}
			for(int i=0; i<row.length; i++){
				require(Double.isFinite(Double.parseDouble(row[i])), "Unexpected nonfinite field "+names[i]+"="+row[i]);
			}
			System.out.println(name+"\t"+row[6]+"\t"+row[7]+"\t"+row[8]+"\t"+row[9]+"\t"+row[24]+"\tPASS");
		}
	}

	/** Checks literal report expectations independently of JVM assertion flags. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}
