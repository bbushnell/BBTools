package test.barcode;

import java.util.Locale;

import barcode.DemuxData;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import structures.ByteBuilder;

/**
 * Finite A48 numeric-domain checks for DemuxData base-5 barcode encoding.
 * Baseline mode records the positive long overflow for a 28-base all-T barcode;
 * fixed mode requires the 28-base lossless control to remain valid and the
 * 28-base overflowing slice to be rejected before a corrupt A48 value is emitted.
 * @author Brian Bushnell, Jean
 */
public final class DemuxDataA48DomainTest{

	/** Requires mode=baseline|fixed; assertions must be enabled. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, DemuxDataA48DomainTest.class, false);
		try{
			String mode=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(!Parser.parseStatic(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			if(!"baseline".equals(mode) && !"fixed".equals(mode)){
				throw new IllegalArgumentException("Set mode=baseline or mode=fixed.");
			}
			boolean assertions=false;
			assert(assertions=true);
			if(!assertions){throw new IllegalStateException("DemuxData A48 domain regression requires -ea.");}

			final boolean baseline=mode.equals("baseline");
			final long maxCode=roundTrip(MAX_27);
			final long longCode=roundTrip(LOSSLESS_28);
			check(longCode==LOSSLESS_28_CODE, "Unexpected 28-base lossless code: "+longCode);
			final long maxLongCode=roundTrip(LONG_MAX_28);
			check(maxLongCode==Long.MAX_VALUE, "Unexpected Long.MAX_VALUE code: "+maxLongCode);
			expectMinOverflow(baseline);
			boolean rejected=false;
			long wrapped=Long.MIN_VALUE;
			try{
				wrapped=DemuxData.encodeACGTN(OVERFLOW_28, 0, OVERFLOW_28.length());
			}catch(final AssertionError e){rejected=true;}
			catch(final IllegalArgumentException e){rejected=true;}
			check(rejected!=baseline, "28-base slice rejection mismatch: mode="+mode+", rejected="+rejected);
			if(baseline){
				check(wrapped==WRAPPED_28_T, "Unexpected wrapped 28-base code: "+wrapped);
				final ByteBuilder decoded=new ByteBuilder();
				DemuxData.decodeACGTN(wrapped, decoded, OVERFLOW_28.length());
				check(!OVERFLOW_28.equals(decoded.toString()), "Wrapped 28-base code must not round-trip silently.");
				System.out.println("wrapped_28T="+wrapped+" decoded="+decoded);
			}
			System.out.println("DEMUXDATA_A48_DOMAIN mode="+mode+" max27="+maxCode+
						" lossless28="+longCode+" maxLong28="+maxLongCode+" rejected28="+rejected);
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** No instances; each process owns the bounded representable/overflow domain check. */
	private DemuxDataA48DomainTest(){}

	/** Round-trips one legal A/C/G/N/T slice through DemuxData's base-5 codec. */
	private static long roundTrip(final String barcode){
		final long code=DemuxData.encodeACGTN(barcode, 0, barcode.length());
		final ByteBuilder decoded=new ByteBuilder();
		DemuxData.decodeACGTN(code, decoded, barcode.length());
		check(barcode.equals(decoded.toString()), "A48 round trip failed: "+barcode+" -> "+code+" -> "+decoded);
		return code;
	}

	/**
	 * Verifies that the immediate unsigned value after Long.MAX_VALUE is rejected.
	 * Baseline code rejects it by assertion after overflow; fixed code rejects it
	 * before overflow by throwing IllegalArgumentException.
	 */
	private static void expectMinOverflow(final boolean baseline){
		Throwable thrown=null;
		try{
			DemuxData.encodeACGTN(UNSIGNED_MIN_28, 0, UNSIGNED_MIN_28.length());
		}catch(final AssertionError e){
			thrown=e;
		}catch(final IllegalArgumentException e){
			thrown=e;
		}
		check(thrown!=null, "Unsigned 2^63 value encoded without rejection.");
		if(baseline){
			check(thrown instanceof AssertionError,
					"Baseline should reject unsigned 2^63 by assertion: "+thrown);
		}else{
			check(thrown instanceof IllegalArgumentException,
					"Fixed mode should reject unsigned 2^63 by IllegalArgumentException: "+thrown);
		}
	}

	/** Always enforce fixture oracles independent of assertion configuration. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}

	private static final String MAX_27="TTTTTTTTTTTTTTTTTTTTTTTTTTT";
	private static final String LOSSLESS_28="A"+MAX_27;
	private static final String OVERFLOW_28=MAX_27+"T";
	private static final String LONG_MAX_28="CCATNNGTACNATTGGTNTNCANCCGCG";
	private static final String UNSIGNED_MIN_28="CCATNNGTACNATTGGTNTNCANCCGCN";
	private static final long LOSSLESS_28_CODE=7450580596923828124L;
	private static final long WRAPPED_28_T=359414837200037392L;
}
