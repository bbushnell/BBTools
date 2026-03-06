package cardinality;

/**
 * Interface for cardinality trackers that support calibration and comparison.
 * Implemented by CardinalityTracker (with default UnsupportedOperationException stubs)
 * and overridden by DDL/DDL2/DDL8/DLL4 subclasses.
 * <p>
 * Allows DDLCalibrationDriver to drive any tracker type via a loglogtype parameter
 * without knowing the concrete class.
 *
 * @author Chloe
 * @date March 2026
 */
public interface Drivable {

	/** Returns all cardinality estimates for calibration. */
	double[] rawEstimates();

	/** Returns number of filled (non-empty) buckets. */
	int filledBuckets();

	/** Returns fraction of filled buckets: filledBuckets/buckets. */
	double occupancy();

	/** Returns lastCardinality cache (-1 if stale). */
	long getLastCardinality();

	/** Sets lastCardinality cache (use -1 to mark stale). */
	void setLastCardinality(long val);

}
