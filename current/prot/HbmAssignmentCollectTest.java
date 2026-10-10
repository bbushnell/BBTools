package prot;

/** Independent feasible/infeasible paired-coverage geometries. @author Keqing */
public final class HbmAssignmentCollectTest {
	public static void main(String[] args){
		HbmAssignmentCollect.checkCoverage(80, 100, 0, 99, 0, 79, 1, .8f, .8f);
		HbmAssignmentCollect.checkCoverage(101, 140, 20, 119, 20, 110, (float)(91/101.0), .91f, (float)(91/101.0));
		reject(()->HbmAssignmentCollect.checkCoverage(10, 100, 0, 99, 0, 99, 1, 1, 1), "integer paired-column");
		reject(()->HbmAssignmentCollect.checkCoverage(100, 100, 0, 99, 0, 70, 1, 1, 1), "intersection");
		reject(()->HbmAssignmentCollect.checkCoverage(80, 140, 40, 139, 0, 79, 1, .8f, .8f), "intersection");
		reject(()->HbmAssignmentCollect.checkCoverage(100, 100, 0, 99, 0, 99, Float.NaN, 1, 1), "coverage gate");
		reject(()->HbmAssignmentCollect.checkCoverage(Integer.MAX_VALUE, Integer.MAX_VALUE, 0, Integer.MAX_VALUE-1, 0, Integer.MAX_VALUE-1, 1, 1, 1), "bounds");
		System.err.println("HBM_COLLECTION_GEOMETRY_TEST_PASS feasible=true integer_count=true span_intersection=true malformed_rejected=true");
	}
	private static void reject(Runnable action, String diagnostic){
		try{action.run();}catch(IllegalArgumentException e){if(!e.getMessage().contains(diagnostic)){throw new AssertionError("Unexpected diagnostic: "+e);} return;}
		throw new AssertionError("Impossible geometry was accepted: "+diagnostic);
	}
}
