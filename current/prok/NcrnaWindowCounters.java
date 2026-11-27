package prok;

import java.util.Arrays;

import shared.KillSwitch;

/** Optional primitive observations of a caller's actual model/mask work.
 * Allocated only by an explicit observer; ordinary serving keeps a null hook.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaWindowCounters implements NcrnaAlignmentWork.QueryObserver{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	NcrnaWindowCounters(final int models){
		assert(models>0) : "Kernel observations require the actual nonempty global model axis";
		final int size=Math.multiplyExact(models, 8);
		requests=new long[size];
		physical=new long[size];
		queryBases=new long[size];
		modelBases=new long[size];
		accepted=new long[size];
		min=new int[size];
		max=new int[size];
		Arrays.fill(min, Integer.MAX_VALUE);
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	void activate(final int model, final int mask){
		assert(model>=0 && model<requests.length/8 && mask>0 && mask<8)
			: KillSwitch.assertDie("Kernel attribution borrows a real global model and the scheduled member mask");
		active=model*8+mask;
	}
	@Override
	public void query(final int modelLength, final int queryLength, final boolean executed){
		assert(active>=0 && modelLength>0 && queryLength>0)
			: KillSwitch.assertDie("Nested kernel widths need an actual active model and nonempty operands");
		requests[active]++;
		if(executed){physical[active]++;}
		queryBases[active]+=queryLength;
		modelBases[active]+=modelLength;
		min[active]=Math.min(min[active], queryLength);
		max[active]=Math.max(max[active], queryLength);
	}
	void accepted(final int model, final int mask){
		assert(model>=0 && model<requests.length/8 && mask>0 && mask<8)
			: KillSwitch.assertDie("Accepted proposals retain the actual final model and source-member mask");
		accepted[model*8+mask]++;
	}
	long sum(final long[] values){
		assert(values!=null && values.length==requests.length) : "Kernel populations share one model/mask axis";
		long result=0;
		for(final long value:values){result+=value;}
		return result;
	}
	void add(final NcrnaWindowCounters other){
		assert(other!=null && other!=this && requests.length==other.requests.length)
			: KillSwitch.assertDie("Joined kernel observations must retain the exact bound model/mask axes");
		for(int index=0; index<requests.length; index++){
			requests[index]+=other.requests[index];
			physical[index]+=other.physical[index];
			queryBases[index]+=other.queryBases[index];
			modelBases[index]+=other.modelBases[index];
			accepted[index]+=other.accepted[index];
			min[index]=Math.min(min[index], other.min[index]);
			max[index]=Math.max(max[index], other.max[index]);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public final long[] requests, physical, queryBases, modelBases, accepted;
	public final int[] min, max;
	private int active=-1;
}
