package test.barcode;

import java.util.Locale;

import barcode.PCRMatrixHDist;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;

/**
 * Checks the three aggregate fraction contracts, including an empty matrix.
 * Baseline mode qualifies the existing NaNs; fixed mode requires explicit failure.
 * Positive-denominator controls use independently specified exact binary fractions.
 * @author Brian Bushnell, Shinobu
 */
public final class MatrixFractionRegressionTest{

	/** Requires mode=baseline|fixed and assertions enabled. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, MatrixFractionRegressionTest.class, false);
		try{
			String mode=null;
			for(final String arg : pp.args){
				final String[] pair=Parse.splitOnFirst(arg, '=');
				final String key=pair[0].toLowerCase(Locale.ROOT), value=pair.length>1 ? pair[1] : null;
				if(key.equals("mode")){mode=value;}
				else if(!Parser.parseStatic(arg, key, value)){throw new IllegalArgumentException("Unknown argument: "+arg);}
			}
			check("baseline".equals(mode) || "fixed".equals(mode), "Set mode=baseline or mode=fixed.");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Matrix fraction checks require -ea.");
			final Matrix matrix=new Matrix();
			int undefined=0;
			for(int metric=0; metric<3; metric++){
				if("baseline".equals(mode)){
					check(Float.isNaN(matrix.fraction(metric)), "Empty baseline metric "+metric+" must expose the recorded NaN.");
					undefined++;
				}else{
					boolean threw=false;
					try{matrix.fraction(metric);}
					catch(final IllegalStateException expected){threw=true;}
					check(threw, "Empty metric "+metric+" must reject an undefined denominator.");
					undefined++;
				}
			}
			checkFractions(matrix, 4, 3, 2, 0.75f, 0.5f, 0.25f);
			checkFractions(matrix, 4, 0, 0, 0, 0, 0);
			checkFractions(matrix, 4, 4, 4, 1, 1, 0);
			System.out.println("MATRIX_FRACTIONS mode="+mode+" empty_checks="+undefined+" finite_checks=9");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Tests ratios without coupling this arithmetic contract to candidate search or worker selection. */
	private static void checkFractions(final Matrix matrix, final long counted, final long assigned,
			final long expected, final float assignedFraction, final float expectedFraction, final float chimericFraction){
		matrix.setTotals(counted, assigned, expected);
		check(matrix.assignedFraction()==assignedFraction, "Assigned fraction differs from literal control.");
		check(matrix.expectedFraction()==expectedFraction, "Expected fraction differs from literal control.");
		check(matrix.chimericFraction()==chimericFraction, "Chimeric fraction differs from literal control.");
	}

	/** Exposes only protected aggregate state; fraction implementations remain production methods. */
	private static final class Matrix extends PCRMatrixHDist{
		Matrix(){super(4, 0, 0, false);}
		void setTotals(final long counted, final long assigned, final long expected){
			check(counted>=assigned && assigned>=expected && expected>=0,
					"Fixture totals must represent expected assignments within all assignments within observations.");
			totalCounted=counted;
			totalAssigned=assigned;
			totalAssignedToExpected=expected;
		}
		float fraction(final int metric){
			switch(metric){
				case 0: return assignedFraction();
				case 1: return expectedFraction();
				case 2: return chimericFraction();
				default: throw new IllegalArgumentException("Unknown fraction metric: "+metric);
			}
		}
	}

	/** Enforces independent test predicates regardless of production assertion settings. */
	private static void check(final boolean condition, final String message){
		if(!condition){throw new AssertionError(message);}
	}
	private MatrixFractionRegressionTest(){}
}
