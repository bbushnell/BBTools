package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

import dna.AminoAcid;
import map.LongIntMap;
import parse.LineParser1;
import prok.NcrnaTextIO.Input;
import shared.KillSwitch;
import structures.IntList;
import structures.LongList;

/** Driver-only overlap clustering; one mutable builder belongs to one worker.
 * Raw intervals use reusable primitive scratch. Only returned windows allocate.
 * This constructs geometry, not seed support: the caller still applies its
 * unchanged distinct-key admission to the final pooled query. Across-member
 * pooling belongs to NcrnaClanScavenger and must always use union bounds.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaOverlapWindowBuilder{

	/** Blanket bounds preserve the ordinary center/padding/K convention. Per-key
	 * bounds instead extend the actual inclusive seventeen-base seed footprint,
	 * exactly as ReportNcrnaSeedPositions.Scratch.proposals does. */
	public ArrayList<int[]> build(final int[] centers, final long[] keys, final int sequenceLength,
		final int padding, final boolean union, final int minimum, final int edgeMinimum, final Caps caps){
		return build(centers, keys, sequenceLength, padding, union, minimum, edgeMinimum, caps, 0);
	}

	/** A positive ceiling splits overlap chains before repeated-key enveloping.
	 * Original single-key intervals are retained even when wider than the ceiling. */
	public ArrayList<int[]> build(final int[] centers, final long[] keys, final int sequenceLength,
		final int padding, final boolean union, final int minimum, final int edgeMinimum, final Caps caps, final int ceiling){
		require(ceiling>=0 && (ceiling==0 || union && caps!=null),
			"Chain ceilings apply only to the explicit per-key union experiment");
		require(centers!=null && sequenceLength>0 && padding>=0 && minimum>0 && edgeMinimum>=0
			&& keys!=null && keys.length==centers.length,
			"Overlap windows require actual seed centers, positive lengths and the identity of every seed occurrence");
		clearGeometry();
		unavailableSeeds=0;
		ensure(centers.length);
		int n=0;
		for(int i=0; i<centers.length; i++){
			final int center=centers[i];
			require(center>=K/2 && (long)center+K/2<sequenceLength
				&& (i==0 || center>centers[i-1]), "Seed centers must be unique, ordered and contain their complete seventeen-base footprint");
			final long start, stop;
			if(caps==null){
				start=(long)center-padding;
				stop=(long)center+padding+K;
			}else{
				final int index=caps.index.get(keys[i]);
				if(index<0){throw new IllegalArgumentException("The configured per-key bounds omit a scanned seed: "+keys[i]);}
				if(caps.left[index]<0){unavailableSeeds++;continue;}
				final long keyStart=(long)center-K/2;
				start=keyStart-caps.left[index];
				stop=keyStart+K-1+caps.right[index];
			}
			final int a=(int)Math.max(0L, start), z=(int)Math.min(sequenceLength-1L, stop);
			assert(a>=0 && z>=a && z<sequenceLength) : "A nonnegative seed margin must retain its physical seed after contig clipping";
			rawStops[n]=z;
			rawKeys[n]=keyIndex(keys[i]);
			order[n]=((long)a<<32)|(n&0xffffffffL);
			n++;
			geometry.rawWindows++;
			geometry.rawBases+=z-(long)a+1;
		}
		Arrays.sort(order, 0, n);
		final ArrayList<int[]> result=new ArrayList<int[]>();
		if(n==0){return result;}
		for(int first=0; first<n; ){
			int end=first+1, lastStop=rawStops[(int)order[first]];
			final int start=(int)(order[first]>>>32);
			int limit=ceiling==0 ? 0 : Math.max(ceiling, lastStop-start+1);
			while(end<n && (int)(order[end]>>>32)<=lastStop){
				final int nextStart=(int)(order[end]>>>32), nextStop=rawStops[(int)order[end]];
				final int stop=Math.max(lastStop, nextStop);
				final int nextLimit=ceiling==0 ? 0 : Math.max(limit, nextStop-nextStart+1);
				if(ceiling>0 && stop-(long)start+1>nextLimit){break;}
				lastStop=stop;
				limit=nextLimit;
				end++;
			}
			combine(first, end, union, sequenceLength, minimum, edgeMinimum, limit, result);
			first=end;
		}
		assert(geometry.conserves()) : "Overlap groups must conserve each usable seed occurrence exactly once, including discarded groups";
		return result;
	}

	/** First union repeated occurrences of each key inside this connected raw
	 * overlap component. Only DISTINCT-key envelopes take part in intersections.
	 * A key repeated in disjoint components cannot bridge the gap by itself. */
	private void combine(final int first, final int end, final boolean union, final int length,
		final int minimum, final int edgeMinimum, final int ceiling, final ArrayList<int[]> output){
		assert(first>=0 && end>first && end<=order.length)
			: "A component owns a nonempty sorted slice of actual raw seed occurrences";
		active.clear();
		for(int i=first; i<end; i++){
			final int raw=(int)order[i], key=rawKeys[raw], a=(int)(order[i]>>>32), z=rawStops[raw];
			if(counts[key]==0){
				active.add(key);
				left[key]=a;
				right[key]=z;
			}else{
				left[key]=Math.min(left[key], a);
				right[key]=Math.max(right[key], z);
			}
			counts[key]++;
		}
		envelopes.clear();
		for(int i=0; i<active.size; i++){
			final int key=active.get(i);
			envelopes.add(((long)left[key]<<32)|(key&0xffffffffL));
		}
		Arrays.sort(envelopes.array, 0, envelopes.size);
		int key=(int)envelopes.get(0), start=left[key], stop=right[key], contributors=counts[key];
		for(int i=1; i<envelopes.size; i++){
			key=(int)envelopes.get(i);
			final int a=left[key], z=right[key];
			if(a<=stop && z>=start){
				start=union ? Math.min(start, a) : Math.max(start, a);
				stop=union ? Math.max(stop, z) : Math.min(stop, z);
				contributors+=counts[key];
			}else{
				emit(start, stop, contributors, length, minimum, edgeMinimum, ceiling, output);
				start=a;
				stop=z;
				contributors=counts[key];
			}
		}
		emit(start, stop, contributors, length, minimum, edgeMinimum, ceiling, output);
		for(int i=0; i<active.size; i++){counts[active.get(i)]=0;}
	}
	/** The fixed scanner seed set bounds this reusable primitive dictionary;
	 * only the first occurrence of a previously unseen key grows its slot state. */
	private int keyIndex(final long word){
		require(word>=0 && word<(1L<<(2*K)), "Geometry must retain a valid actual forward seventeen-mer identity");
		int index=keySlots.get(word);
		if(index>=0){return index;}
		index=keySlots.size();
		if(index==counts.length){
			final int size=Math.multiplyExact(counts.length, 2);
			left=Arrays.copyOf(left, size);
			right=Arrays.copyOf(right, size);
			counts=Arrays.copyOf(counts, size);
		}
		keySlots.put(word, index);
		return index;
	}

	/** Preserve the existing inclusive-length floor and physical-edge exception.
	 * The historical stop-start rather than stop-start+1 test is deliberate here:
	 * changing that separate caller policy would confound this experiment. */
	private void emit(final int start, final int stop, final int contributors, final int length,
		final int minimum, final int edgeMinimum, final int ceiling, final ArrayList<int[]> output){
		assert(start>=0 && stop>=start && stop<length && contributors>0)
			: "Every emitted or discarded overlap group must own a nonempty physical interval and contributing seeds";
		final long width=stop-(long)start+1;
		assert(ceiling==0 || width<=ceiling) : "Repeated-key envelopes cannot undo the bounded raw-component split";
		final boolean retained=stop-start>=minimum
			|| edgeMinimum>0 && (start==0 || stop==length-1) && width>=edgeMinimum;
		if(!retained){
			geometry.droppedWindows++;
			geometry.droppedContributors+=contributors;
		}else{
			if(contributors==1){
				geometry.separateWindows++;
				geometry.separateBases+=width;
			}else{
				// Existing counter storage; the new experiment must label this as
				// combined groups, not claim union groups were intersections.
				geometry.intersectedWindows++;
				geometry.intersectedBases+=width;
				geometry.intersectedContributors+=contributors;
			}
			output.add(ceiling==0 ? new int[]{start, stop, NcrnaScavenger.WINDOW_COLLAPSE}
				: new int[]{start, stop, NcrnaScavenger.WINDOW_COLLAPSE, ceiling});
		}
	}
	private void ensure(final int n){
		assert(n>=0) : "Primitive overlap scratch is sized from an actual seed-position array";
		if(order.length<n){
			final int size=Math.max(n, Math.max(16, order.length*2));
			order=KillSwitch.allocLong1D(size);
			rawStops=KillSwitch.allocInt1D(size);
			rawKeys=KillSwitch.allocInt1D(size);
		}
	}
	private void clearGeometry(){
		geometry.rawWindows=0;
		geometry.rawBases=0;
		geometry.separateWindows=0;
		geometry.separateBases=0;
		geometry.intersectedWindows=0;
		geometry.intersectedBases=0;
		geometry.intersectedContributors=0;
		geometry.droppedWindows=0;
		geometry.droppedContributors=0;
	}

	/** Owns mutable proposals for one physical read/strand/pass. Across members
	 * only union is valid: a seed does not identify the eventual winning member.
	 * No seed-count or model-eligibility decision is made here. */
	public static ArrayList<int[]> pool(final ArrayList<int[]> proposals){
		require(proposals!=null && (proposals.isEmpty() || proposals.get(0)!=null),
			"Pool only one actual physical read/strand/pass schedule");
		final int fields=proposals.isEmpty() ? 4 : proposals.get(0).length;
		for(final int[] window:proposals){
			require(window!=null && (fields==4 || fields==5) && window.length==fields && window[0]>=0 && window[1]>=window[0]
				&& window[2]==NcrnaScavenger.WINDOW_COLLAPSE && window[3]>0 && window[3]<8,
				"Pooled proposals require physical bounds, ordinary source and a nonempty three-member mask");
			require(fields==4 || window[4]>=window[1]-(long)window[0]+1,
				"A protected proposal must retain its member ceiling or its wider original interval");
		}
		proposals.sort(BY_SPAN);
		final ArrayList<int[]> result=new ArrayList<int[]>();
		for(final int[] window:proposals){
			final int[] last=result.isEmpty() ? null : result.get(result.size()-1);
			final int limit=fields==5 && last!=null ? Math.max(last[4], window[4]) : 0;
			if(last!=null && window[0]<=last[1]
				&& (fields==4 || Math.max(last[1], window[1])-(long)last[0]+1<=limit)){
				last[1]=Math.max(last[1], window[1]);
				last[3]|=window[3];
				if(fields==5){last[4]=limit;}
			}else{result.add(window);}
		}
		return result;
	}

	/** Immutable, occurrence-derived margins. A known unobserved key has -1 at
	 * both ends; an entirely missing key is a configuration error. Keys stay on
	 * the explicit forward axis; they are not canonicalized or reverse-complemented. */
	public static final class Caps{
		/** Read one complete member's fixed-axis p99+10 table once, before workers
		 * start. Legacy p99.9 tables are intentionally not interchangeable. */
		public static Caps load(final String path, final String member, final int expectedKeys){
			return load(path, member, expectedKeys, 99);
		}
		/** An explicit percentile must match every observed table row; omission
		 * above keeps the historical p99 contract, including its diagnostics. */
		public static Caps load(final String path, final String member, final int expectedKeys, final int percentile){
			require(path!=null && member!=null && expectedKeys>0, "Per-key bounds need an explicit member, table and fixed-axis size");
			final String rule=capRule(percentile);
			final LongList words=new LongList(expectedKeys);
			final IntList left=new IntList(expectedKeys), right=new IntList(expectedKeys);
			final LineParser1 p=new LineParser1('\t');
			long line=1;
			try(Input input=new Input(path)){
				RrnaResourceIO.header(input.reader, CAP_HEADER, path);
				for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
					line++;
					p.set(row);
					require(p.terms()==8 && p.termEquals(member, 0), "Each seed-bound row must retain its exact member and eight-field schema");
					final long word=word(p);
					final long occurrences=p.parseLong(2), references=p.parseLong(3);
					require(occurrences>=0 && references>=0 && references<=occurrences && (occurrences==0)==(references==0),
						"Each reference contributes at least one occurrence; zero-occurrence keys have zero supporting references");
					if(occurrences==0){
						require(p.termEquals("NA", 4) && p.termEquals("NA", 5) && p.termEquals("NA", 6)
							&& p.termEquals("NO_OCCURRENCES", 7), "Unobserved keys retain explicit unavailable bounds");
						left.add(-1);
						right.add(-1);
					}else{
						final int a=p.parseInt(4), z=p.parseInt(5);
						require(a>=0 && z>=0 && p.parseLong(6)==K+(long)a+z && p.termEquals(rule, 7),
							"The per-key caller requires p"+percentile+"+10 bounds capped by the full member maximum-minus17, not legacy p99.9 margins");
						left.add(a);
						right.add(z);
					}
					words.add(word);
				}
			}catch(RuntimeException failure){throw new IllegalArgumentException("Invalid per-key bounds at line "+line+" in "+path, failure);}
			require(words.size==expectedKeys, "Every fixed seed key must retain one available or explicitly unavailable bounds row");
			return new Caps(Arrays.copyOf(words.array, words.size), Arrays.copyOf(left.array, left.size), Arrays.copyOf(right.array, right.size));
		}
		private static long word(final LineParser1 p){
			p.setBounds(1);
			require(p.b()-p.a()==K, "Per-key bounds retain exact forward seventeen-mers");
			long word=0;
			for(int i=p.a(); i<p.b(); i++){
				final int base=p.line()[i]&255;
				final int lower=base|32;
				final int value=base<AminoAcid.baseToNumber.length ? AminoAcid.baseToNumber[base] : -1;
				require(value>=0 && (lower=='a' || lower=='c' || lower=='g' || lower=='t'),
					"A fixed forward seed word contains only A/C/G/T bases");
				word=(word<<2)|value;
			}
			return word;
		}
		public Caps(final long[] keys, final int[] left_, final int[] right_){
			require(keys!=null && left_!=null && right_!=null && keys.length>0
				&& keys.length==left_.length && keys.length==right_.length,
				"Per-key margins require one left/right row on the same nonempty fixed forward seed axis");
			left=left_.clone();
			right=right_.clone();
			index=new LongIntMap(keys.length*2);
			for(int i=0; i<keys.length; i++){
				require(keys[i]>=0 && keys[i]<(1L<<(2*K)) && index.get(keys[i])<0,
					"The fixed forward seventeen-mer axis must contain unique valid words");
				require(left[i]>=-1 && right[i]>=-1 && (left[i]<0)==(right[i]<0),
					"Unavailable margins use -1 at both ends; observed margins must be nonnegative");
				index.put(keys[i], i);
			}
		}
		/** A row count cannot establish the loaded forward seed identity. Check
		 * the complete actual resource set once, before constructing workers. */
		public void validateKeys(final long[] keys){
			require(keys!=null && keys.length==left.length,
				"Per-key margins must use the complete actual loaded forward seed axis");
			for(final long key:keys){
				require(index.get(key)>=0, "Per-key margins omit a word from the actual loaded forward seed axis: "+key);
			}
		}
		/** Diagnostic lookup: -2 means absent from the fixed axis, -1 unavailable. */
		public int leftMargin(final long key){
			final int slot=index.get(key);
			return slot<0 ? -2 : left[slot];
		}
		/** Diagnostic lookup on the same immutable forward axis as leftMargin. */
		public int rightMargin(final long key){
			final int slot=index.get(key);
			return slot<0 ? -2 : right[slot];
		}
		private final LongIntMap index;
		private final int[] left, right;
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	private long[] order=new long[0];
	private int[] rawStops=new int[0], rawKeys=new int[0];
	private int[] left=new int[16], right=new int[16], counts=new int[16];
	private final LongIntMap keySlots=new LongIntMap(16);
	private final IntList active=new IntList();
	private final LongList envelopes=new LongList();
	public final NcrnaScavenger.WindowGeometry geometry=new NcrnaScavenger.WindowGeometry();
	public long unavailableSeeds;
	private static final int K=17;
	public static final String CAP_HEADER="member\tforward_17mer\toccurrences\tdistinct_references\tleft_cap_bp\tright_cap_bp\tper_hit_window_bp\testimator";
	public static final String CAP_RULE="P99_PLUS10_MEMBER_MAX_MINUS17";
	/** Shared producer/consumer label; only the four measured percentile rules
	 * are supported, and a mislabeled or mixed table fails before callers start. */
	public static String capRule(final int percentile){
		require(percentile==90 || percentile==95 || percentile==98 || percentile==99,
			"Bound percentile must be one of the measured 90,95,98,99 rules");
		return percentile==99 ? CAP_RULE : "P"+percentile+"_PLUS10_MEMBER_MAX_MINUS17";
	}
	private static final Comparator<int[]> BY_SPAN=(a, b)->{
		final int start=Integer.compare(a[0], b[0]);
		return start==0 ? Integer.compare(a[1], b[1]) : start;
	};
}
