package prok;

import java.util.ArrayList;

import map.LongHashSet;
import map.LongIntMap;
import map.IntObjectMap;
import map.ObjectMap;
import prok.NcrnaReferenceCalls.Rule;
import structures.IntList;

/** Exact reported-call losses outside the explicitly included target family.
 * Register one contig's immutable candidate reports, then count selected native
 * candidate IDs. Equivalent alternatives share a coordinate/family group, but
 * selected duplicates retain their multiplicity. A model/ID change at identical
 * coordinates is retained; an endpoint change is an exact-call loss plus gain,
 * not a claim of biological gene loss. No overlap matching or net-count proxy.
 * The caller checks source/snapshot ownership before passing selected IDs.
 * One instance owns one contig and one named baseline; trials reuse its maps,
 * groups and counters without allocating per selected call or trial.
 * @author Brian, Raiden
 */
public final class NcrnaRnaCallChanges{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaRnaCallChanges(final NcrnaReferenceCalls projection_){
		require(projection_!=null, "Other-RNA accounting needs the same explicit target-family projection as RNA grading");
		projection=projection_;
	}

	/** Register potential zero-based, inclusive genomic reports before begin().
	 * CDS and explicitly included target RNAs are ignored. All other structural
	 * and generic RNAs retain their exact type, family, strand and coordinates.
	 */
	public void candidate(final long id, final int type, final int strand, final int start,
			final int stop, final String model, final String family){
		require(!sealed && id>=0 && type>=ProkObject.CDS && type<=ProkObject.RNA
			&& strand>=0 && strand<2 && start>=0 && stop>=start,
			"Candidate registration precedes selection and requires valid native ID/type/strand/report bounds");
		if(type==ProkObject.CDS){return;}
		final String name;
		if(type==ProkObject.RNA){
			final Rule rule=projection.modelRule(model);
			require(family==null || rule.displayLabel.equals(NcrnaFamilyRegistry.displayName(family)),
				"Other-RNA candidate family must agree with its explicit model projection");
			if(rule.include){return;}
			name=rule.family;
		}else{name=family==null ? ProkObject.typeStrings[type] : NcrnaFamilyRegistry.get(family).id;}
		if(!seen.add(id)){
			throw new IllegalArgumentException("An other-RNA candidate ID must have one immutable report: "+id);
		}
		LongIntMap[] partitions=byFamily.get(name);
		if(partitions==null){
			partitions=new LongIntMap[2*(ProkObject.RNA+1)];
			byFamily.put(name, partitions);
		}
		final int partition=2*type+strand;
		LongIntMap coordinates=partitions[partition];
		if(coordinates==null){
			coordinates=new LongIntMap(8);
			partitions[partition]=coordinates;
		}
		final long key=((long)start<<32)|(stop&0xffffffffL);
		int group=coordinates.get(key);
		if(group<0){
			group=groups.size();
			coordinates.put(key, group);
			groups.add(new Group(type, strand, start, stop));
		}
		candidateGroups.put(id, group);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Selections          ----------------*/
	/*--------------------------------------------------------------*/

	/** Start one selection; candidate registration is permanently closed. */
	public void begin(){
		require(!active, "Finish the previous RNA selection before beginning another");
		sealed=true;
		active=true;
		ready=false;
		selected=0;
		seen.clear();
		for(int i=0; i<groups.size(); i++){groups.get(i).current=0;}
	}
	/** Count a selected ID from the enclosing verified snapshot. IDs outside the
	 * registered other-RNA pool are CDS or included-target calls and are ignored.
	 * Two distinct selected IDs at identical coordinates count twice.
	 */
	public void select(final long id){
		require(active && id>=0, "RNA selection IDs need an active selection and a validated native candidate identity");
		final int index=candidateGroups.get(id);
		if(index<0){return;}
		if(!seen.add(id)){
			throw new IllegalArgumentException("A completed native path cannot repeat an other-RNA candidate ID: "+id);
		}
		groups.get(index).current++;
		selected++;
	}
	/** Checked native-report route. Candidate IDs are still source/snapshot-bound
	 * by the caller; structural RNAs also validate their frozen report geometry.
	 */
	public void select(final long id, final int type, final int strand, final int start, final int stop){
		require(active && id>=0, "Report selection requires an active trial and a validated snapshot ID");
		final int index=candidateGroups.get(id);
		if(index>=0){
			final Group group=groups.get(index);
			require(group.type==type && group.strand==strand && group.start==start && group.stop==stop,
				"Selected other-RNA reporting geometry must match its immutable candidate projection");
		}
		select(id);
	}
	/** Finish counts; freezeBaseline(name) explicitly chooses the comparison arm. */
	public void finish(){
		require(active, "RNA selection must begin before its counts can finish");
		active=false;
		ready=true;
	}
	/** Bind the current completed selection exactly once, including an empty one. */
	public void freezeBaseline(final String name){
		require(ready && baseline==null && name!=null && !name.isEmpty(),
			"Freeze one completed, explicitly named other-RNA baseline before comparing trials");
		baseline=name;
		baselineCalls=selected;
		for(int i=0; i<groups.size(); i++){
			final Group group=groups.get(i);
			group.baseline=group.current;
		}
	}
	public String baseline(){check(); return baseline;}
	public int baselineCalls(){check(); return baselineCalls;}
	public int selectedCalls(){check(); return selected;}
	/** Every lost coordinate/family copy is penalized, regardless of gains elsewhere. */
	public int lost(){
		check();
		int lost=0;
		for(int i=0; i<groups.size(); i++){
			final Group group=groups.get(i);
			lost+=Math.max(0, group.baseline-group.current);
		}
		return lost;
	}
	public int gained(){
		check();
		int gained=0;
		for(int i=0; i<groups.size(); i++){
			final Group group=groups.get(i);
			gained+=Math.max(0, group.current-group.baseline);
		}
		return gained;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Helpers           ----------------*/
	/*--------------------------------------------------------------*/

	private void check(){require(ready && baseline!=null, "Other-RNA changes require a complete trial and its explicit frozen baseline");}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/** Exact multiset comparison of two saved-GFF arms for ONE physical contig.
	 * The caller assigns a shared nonnegative partition to each family/type/strand
	 * and omits the target families. Coordinates use the same inclusive convention
	 * in both arms. Models and row IDs are deliberately not part of call identity.
	 * Primitive groups retain multiplicity without allocating an object per row.
	 * One thread owns an instance; finish once before reading its totals.
	 */
	public static final class Reported{
		/** Add one reported copy to arm0(before) or arm1(after), before finish(). */
		public void add(final int partition, final int start, final int stop, final int arm){
			require(!finished && partition>=0 && start>=0 && stop>=start && arm>=0 && arm<2,
				"Saved RNA groups need a valid family/type/strand partition and equal-coordinate convention before finish");
			LongIntMap coordinates=partitions.get(partition);
			if(coordinates==null){
				coordinates=new LongIntMap(8);
				partitions.put(partition, coordinates);
			}
			final long key=((long)start<<32)|(stop&0xffffffffL);
			int group=coordinates.get(key);
			if(group<0){
				group=counts[0].size;
				coordinates.put(key, group);
				counts[0].add(0);
				counts[1].add(0);
			}
			final IntList selected=counts[arm];
			selected.set(group, Math.addExact(selected.get(group), 1));
		}
		/** Zero calls are a measured empty pair, not missing data. */
		public void finish(){
			require(!finished, "Finish each saved RNA pair once so totals cannot be accumulated twice");
			for(int i=0; i<counts[0].size; i++){
				final int old=counts[0].get(i), next=counts[1].get(i);
				before=Math.addExact(before, old);
				after=Math.addExact(after, next);
				lost=Math.addExact(lost, Math.max(0, old-next));
				gained=Math.addExact(gained, Math.max(0, next-old));
			}
			assert(before-lost+gained==after) : "Exact saved RNA multiset groups must conserve before-loss+gain=after";
			finished=true;
		}
		public long before(){checkFinished(); return before;}
		public long after(){checkFinished(); return after;}
		public long lost(){checkFinished(); return lost;}
		public long gained(){checkFinished(); return gained;}
		private void checkFinished(){require(finished, "Read saved RNA counts only after both complete arms finish");}
		private final IntObjectMap<LongIntMap> partitions=new IntObjectMap<LongIntMap>(8, LongIntMap.class);
		private final IntList[] counts={new IntList(16), new IntList(16)};
		private long before, after, lost, gained;
		private boolean finished;
	}
	private static final class Group{
		Group(final int type_, final int strand_, final int start_, final int stop_){
			type=type_;
			strand=strand_;
			start=start_;
			stop=stop_;
		}
		final int type, strand, start, stop;
		int baseline, current;
	}
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	private final NcrnaReferenceCalls projection;
	private final LongIntMap candidateGroups=new LongIntMap(16);
	private final ObjectMap<String,LongIntMap[]> byFamily=new ObjectMap<String,LongIntMap[]>(String.class, LongIntMap[].class);
	private final ArrayList<Group> groups=new ArrayList<Group>();
	/** Registration inserts the full candidate pool, sizing storage before trials. */
	private final LongHashSet seen=new LongHashSet(16);
	private boolean sealed, active, ready;
	private int selected, baselineCalls;
	private String baseline;
}
