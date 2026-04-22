package cardinality;

/**
 * Validates multi-way merge accuracy for UDLL6 and ErtlULL.
 * Splits N sequential elements across partitions, merges back,
 * and compares against a single-stream reference estimator.
 * Also tests the master+worker merge pattern used by LogLogWrapper.
 *
 * @author Brian Bushnell
 */
public class MergeTest {

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		final long N=10000000;
		final int splits=16;

		// UDLL6 multi-way merge
		final UltraDynamicLogLog6 ref=new UltraDynamicLogLog6(2048, 31, -1, 0);
		final UltraDynamicLogLog6[] parts=new UltraDynamicLogLog6[splits];
		for(int s=0; s<splits; s++){parts[s]=new UltraDynamicLogLog6(2048, 31, -1, 0);}
		for(long i=0; i<N; i++){
			ref.add(i);
			parts[(int)(i%splits)].add(i);
		}
		for(int s=1; s<splits; s++){parts[0].add(parts[s]);}
		CardinalityTracker.clampToAdded=false;
		System.out.printf("UDLL6 16-way merge: ref=%d merged=%d ratio=%.4f%n",
			ref.cardinality(), parts[0].cardinality(), (double)parts[0].cardinality()/ref.cardinality());

		// ErtlULL multi-way merge
		final ErtlULL eref=new ErtlULL(2048, 31, -1, 0);
		final ErtlULL[] eparts=new ErtlULL[splits];
		for(int s=0; s<splits; s++){eparts[s]=new ErtlULL(2048, 31, -1, 0);}
		for(long i=0; i<N; i++){
			eref.add(i);
			eparts[(int)(i%splits)].add(i);
		}
		for(int s=1; s<splits; s++){eparts[0].add(eparts[s]);}
		System.out.printf("ErtlULL 16-way merge: ref=%d merged=%d ratio=%.4f%n",
			eref.cardinality(), eparts[0].cardinality(), (double)eparts[0].cardinality()/eref.cardinality());

		// Master+worker pattern (like LogLogWrapper)
		final UltraDynamicLogLog6 master=new UltraDynamicLogLog6(2048, 31, -1, 0);
		final UltraDynamicLogLog6 worker=new UltraDynamicLogLog6(2048, 31, -1, 0);
		for(long i=0; i<N; i++){worker.add(i);}
		master.add(worker);
		System.out.printf("UDLL6 master+worker: ref=%d merged=%d ratio=%.4f (master mz=%d worker mz=%d)%n",
			ref.cardinality(), master.cardinality(), (double)master.cardinality()/ref.cardinality(),
			master.getMinZeros(), worker.getMinZeros());
	}

}
