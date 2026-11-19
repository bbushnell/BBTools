package prok;

import java.util.Arrays;

/** Worker-owned nucleotide assignment: DDL order changes the winner only in
 * the explicitly selected top-window or patience policies. No shortlist gate
 * can create a new family before every existing model has failed.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaClusterAssignment{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaClusterAssignment(final Policy policy_, final int patience_, final int bins_){
		if(policy_==null || patience_<0 || bins_<1){
			throw new IllegalArgumentException("Assignment needs an explicit policy, nonnegative patience and positive sketch width");
		}
		policy=policy_;
		patience=patience_;
		bins=bins_;
		histogram=new int[bins+1];
	}

	/*--------------------------------------------------------------*/
	/*----------------          Assignment          ----------------*/
	/*--------------------------------------------------------------*/

	/** Counts and model slots share stable IDs. Null models are inactive; zero
	 * counts remain eligible. Alignment state belongs to this caller's worker. */
	public int bestFamily(final byte[] query, final byte[][] models, final int[] counts,
			final float minimumIdentity, final Alignment alignment){
		if(query==null || query.length==0 || models==null || counts==null || counts.length<models.length
				|| alignment==null || !Float.isFinite(minimumIdentity) || minimumIdentity<0 || minimumIdentity>1){
			throw new IllegalArgumentException("Assignment needs a real query, complete stable model/count vectors and a finite identity gate");
		}
		queries++;
		final int size=rank(models, counts);
		int best=-1, stale=0;
		float bestIdentity=-1;
		for(int rank=0; rank<size; rank++){
			if(policy==Policy.TOP50 && rank==50 && best>=0){break;}
			final int family=order[rank];
			final float identity=alignment.identity(query, models[family]);
			alignments++;
			if(!Float.isFinite(identity) || identity<0 || identity>1){
				throw new IllegalStateException("Nucleotide alignment identity must be finite in [0,1]: family="+family+", identity="+identity);
			}
			final boolean hadWinner=best>=0;
			final boolean passed=identity>=minimumIdentity;
			final boolean improved=passed && identity>bestIdentity;
			if(passed){
				passing++;
				if(improved || identity==bestIdentity && family<best){
					best=family;
					bestIdentity=identity;
				}
			}
			if(improved){stale=0;}
			else if(hadWinner){stale++;}
			if(policy==Policy.PATIENCE && best>=0 && (patience==0 || stale>=patience)){
				if(rank+1<size){earlyStops++;}
				break;
			}
		}
		if(best<0){unassigned++;}
		assert(best<0 || bestIdentity>=minimumIdentity)
			: "Only an actually passing nucleotide alignment may own a member; DDL agreement never grants acceptance";
		return best;
	}

	/** Counting sort keeps descending hit count and ascending stable ID, including
	 * every zero-hit row. The exhaustive arm never pays for DDL ordering. */
	private int rank(final byte[][] models, final int[] counts){
		if(order.length<models.length){order=new int[models.length];}
		int size=0;
		if(policy==Policy.ALL){
			for(int id=0; id<models.length; id++){
				if(models[id]!=null){order[size++]=id;}
			}
			return size;
		}
		Arrays.fill(histogram, 0);
		for(int id=0; id<models.length; id++){
			if(models[id]==null){continue;}
			if(counts[id]<0 || counts[id]>bins){
				throw new IllegalArgumentException("One DDL bin supplies at most one matching vote: family="+id+", hits="+counts[id]);
			}
			histogram[counts[id]]++;
		}
		for(int hits=bins; hits>=0; hits--){
			final int count=histogram[hits];
			histogram[hits]=size;
			size+=count;
		}
		for(int id=0; id<models.length; id++){
			if(models[id]!=null){order[histogram[counts[id]]++]=id;}
		}
		assert(size<=models.length) : "Ranking must conserve the live model roster, including empty sketches";
		return size;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Interfaces          ----------------*/
	/*--------------------------------------------------------------*/

	/** Identity only; the consensus builder retains its separate necessary traces. */
	public interface Alignment{
		float identity(byte[] query, byte[] reference);
	}

	public enum Policy{
		TOP50, PATIENCE, ALL;

		public static Policy parse(final String value){
			if("top50".equalsIgnoreCase(value)){return TOP50;}
			if("patience".equalsIgnoreCase(value)){return PATIENCE;}
			if("all".equalsIgnoreCase(value)){return ALL;}
			throw new IllegalArgumentException("clusterpolicy must be top50, patience or all: "+value);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public long queries, alignments, passing, unassigned, earlyStops;
	private final Policy policy;
	private final int patience, bins;
	private final int[] histogram;
	private int[] order=new int[0];
}
