package prok;

import java.io.PrintStream;
import java.util.ArrayList;

/** One explicit mapping/verifier policy for independently named generic ncRNA members.
 * Reuses the sweep tuple reader; it neither enables callers nor selects defaults.
 * Immutable after loading, applied during resource construction before workers.
 * Table rows take precedence over scalar shortlist flags for the named members;
 * mapping K, seeds, resources, identity, endpoints and path scores are independent.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCallerPolicy{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Loads exactly one named policy; every row names a canonical caller member. */
	public static NcrnaCallerPolicy load(final String path){
		if(path==null || path.isEmpty()){
			throw new IllegalArgumentException("ncrnapolicy requires a joint-policy TSV path");
		}
		final NcrnaJointPolicy[] policies=NcrnaJointPolicy.load(path);
		if(policies.length!=1){
			throw new IllegalArgumentException("Ordinary calling requires exactly one policy, not a sweep: "+path);
		}
		return new NcrnaCallerPolicy(policies[0]);
	}

	private NcrnaCallerPolicy(final NcrnaJointPolicy policy_){
		policy=policy_;
		for(int i=0; i<policy.size(); i++){
			final String name=policy.member(i).name;
			final NcrnaFamilyRegistry.Entry entry=NcrnaFamilyRegistry.get(name);
			if(!entry.id.equals(name)){
				throw new IllegalArgumentException("Policy rows require canonical caller IDs, not aliases or group selectors: "+name);
			}
			if(entry.boundary==null || entry.boundary.trna){
				throw new IllegalArgumentException("Policy rows require a generic shared-boundary caller, not a specialized rRNA/tRNA adapter: "+name);
			}
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Applies only the named member's complete tuple to caller-owned load options. */
	public void apply(final String family, final NcrnaFamilyLoader.Options options){
		assert(family!=null && options!=null)
			: "CallGenes resolves a named family and owns its options before loading shared resources";
		for(int i=0; i<policy.size(); i++){
			final NcrnaJointPolicy.Member member=policy.member(i);
			if(!member.name.equals(family)){continue;}
			options.topN=member.topN;
			options.adaptive=member.adaptive;
			options.floor=member.floor;
			options.minHits=member.minHits;
			options.topFrac=member.topFrac;
			options.qFrac=member.qFrac;
			options.modelPatience=member.patience;
			return;
		}
	}

	/** Every selected row must reach exactly one loaded family with its full tuple.
	 * Runs after all caller/resource flags resolve, before output or worker creation.
	 */
	public void validateLoaded(final ArrayList<NcrnaFamily> families){
		assert(families!=null) : "The ordinary caller publishes its complete family roster before validating explicit policies";
		for(int i=0; i<policy.size(); i++){
			final NcrnaJointPolicy.Member member=policy.member(i);
			NcrnaFamily found=null;
			for(int j=0; j<families.size(); j++){
				final NcrnaFamily family=families.get(j);
				if(!member.name.equals(family.name)){continue;}
				if(found!=null){throw new IllegalArgumentException("Repeated loaded policy member: "+member.name);}
				found=family;
			}
			if(found==null){throw new IllegalArgumentException("Policy member is not enabled/loaded: "+member.name);}
			if(found.indexTopN!=member.topN || found.adaptive!=member.adaptive
				|| found.adaptFloor!=member.floor || found.fixedMinHits!=member.minHits
				|| found.adaptTopFrac!=member.topFrac || found.adaptQFrac!=member.qFrac
				|| found.modelPatience!=member.patience){
				throw new IllegalStateException("Loaded family did not consume its complete explicit policy: "+member.name);
			}
		}
	}

	/** Reports the complete explicit tuples after validateLoaded succeeds; borrows the stream. */
	public void report(final PrintStream out){
		assert(out!=null) : "CallGenes retains ownership of its diagnostic stream";
		for(int i=0; i<policy.size(); i++){
			final NcrnaJointPolicy.Member member=policy.member(i);
			out.println("Explicit ncRNA policy: family="+member.name+" mode="+(member.adaptive ? "adaptive" : "fixed")
				+" topN="+member.topN+" floor="+member.floor+" minHits="+member.minHits
				+" topFrac="+member.topFrac+" qFrac="+member.qFrac+" patience="+member.patience);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final NcrnaJointPolicy policy;
}
