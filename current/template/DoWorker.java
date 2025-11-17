package template;

/**
 * Simple worker interface defining a basic work contract.
 * Implementations must provide the specific work logic through doWork().
 * @author Brian Bushnell
 */
public interface DoWorker {
	
	/** Executes the work assigned to this worker implementation */
	public void doWork();
	
}
