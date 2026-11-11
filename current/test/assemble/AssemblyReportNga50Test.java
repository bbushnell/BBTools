package test.assemble;

import fileIO.ByteFile;
import parse.LineParser1;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/** Checks native reporter summary means and medians against independent literal
 * double-bit oracles. Never recomputes a mean or calls the production mean helper.
 * An incorrect number is printed before the fixture fails. Other report columns
 * and per-genome output bytes are checked by the paired shell recipe.
 * @author Brian Bushnell, Sucrose
 */
public final class AssemblyReportNga50Test{

	/** Requires in=summary.tsv, driver=samek|live and case=large|zero|normal|missing. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, AssemblyReportNga50Test.class, false);
		try{
			final Parser parser=new Parser();
			String driver=null, kind=null;
			for(String arg : pp.args){
				final String[] split=arg.split("=", 2);
				final String key=split[0].toLowerCase(java.util.Locale.ROOT);
				final String value=split.length==2 ? split[1] : null;
				if(key.equals("driver")){driver=value;}
				else if(key.equals("case")){kind=value;}
				else if(!parser.parse(arg, key, value)){throw new IllegalArgumentException("Unknown NGA50 fixture argument: "+arg);}
			}
			if(parser.in1==null || !("samek".equals(driver) || "live".equals(driver))){
				throw new IllegalArgumentException("NGA50 fixture requires in=summary.tsv and driver=samek|live.");
			}
			check(parser.in1, driver, kind);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Reads actual output numbers using Java's exponent-aware parser, independently
	 * of the reporter's decimal formatter. Literal powers of two are exact oracles. */
	private static void check(final String path, final String driver, final String kind){
		final long expected;
		if("large".equals(kind)){expected=0x7fe0000000000000L;}//(H+H)/2=H, H=2^1023
		else if("zero".equals(kind)){expected=0x7fd0000000000000L;}//(H+0)/2=H/2
		else if("normal".equals(kind)){expected=0x4028000000000000L;}//(10+14)/2=12
		else if("missing".equals(kind)){expected=0x4014000000000000L;}//(10+0)/2=5
		else{throw new IllegalArgumentException("Unknown NGA50 fixture case: "+kind);}
		final String[] arms="samek".equals(driver) ? SAMEK_ARMS : LIVE_ARMS;
		final int width="samek".equals(driver) ? 20 : 15;
		final ByteFile input=ByteFile.makeByteFile(path, true);
		final LineParser1 lp=new LineParser1('\t');
		final ByteBuilder report=new ByteBuilder();
		report.append("arm\tfield\texpected_bits\tobserved_bits\tpassed\n");
		int rows=0, failed=0;
		try{
			final byte[] header=input.nextLine();
			if(header==null){throw new IllegalArgumentException("Missing report summary header.");}
			lp.set(header);
			if(lp.terms()!=width || !lp.termEquals("arm", 0) || !lp.termEquals("nga50_zero_mean", 7) ||
					!lp.termEquals("nga50_zero_median", 8)){
				throw new IllegalArgumentException("Wrong report summary schema.");
			}
			for(byte[] line=input.nextLine(); line!=null; line=input.nextLine()){
				lp.set(line);
				if(rows>=arms.length || lp.terms()!=width || !lp.termEquals(arms[rows], 0)){
					throw new IllegalArgumentException("Wrong summary row count, width or arm order: row="+rows);
				}
				for(int column=7; column<=8; column++){
					final double value=Double.parseDouble(lp.parseString(column));
					final long observed=Double.doubleToRawLongBits(value);
					final boolean passed=Double.isFinite(value) && value>=0 && observed==expected;
					if(!passed){failed++;}
					report.append(arms[rows]).tab().append(column==7 ? "mean" : "median").tab()
						.append(Long.toHexString(expected)).tab().append(Long.toHexString(observed)).tab().append(passed).nl();
				}
				rows++;
			}
		}finally{if(input.close()){throw new IllegalStateException("NGA50 fixture summary read failed.");}}
		if(rows!=arms.length){throw new IllegalArgumentException("Missing summary arms: expected="+arms.length+" observed="+rows);}
		report.append("ASSEMBLY_REPORT_NGA50 driver=").append(driver).append(" case=").append(kind)
			.append(" checked=").append(rows*2).append(" failed=").append(failed).nl();
		System.out.print(report.toString());
		System.out.flush();
		if(failed!=0){throw new AssertionError("NGA50 summary numbers differ from literal oracles: "+failed);}
	}

	/** Frozen native default order, checked independently of the report header. */
	private static final String[] SAMEK_ARMS={"baseline", "oracle", "nn_three", "nn_four", "ordinary"};
	private static final String[] LIVE_ARMS={"native", "live"};
}
