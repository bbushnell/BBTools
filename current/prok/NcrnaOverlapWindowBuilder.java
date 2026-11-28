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
		return build(centers, keys, sequenceLength, padding, union, minimum, edgeMinimum, caps, ceiling, false);
	}
	/** New generic files retain unmeasured keys with legacy padding capped at their declared family ceiling. */
	public ArrayList<int[]> build(final int[] centers, final long[] keys, final int sequenceLength,
		final int padding, final boolean union, final int minimum, final int edgeMinimum, final Caps caps,
		final int ceiling, final boolean retainUnmeasured){
		require(ceiling>=0 && (ceiling==0 || caps!=null),
			"Protected overlap chains require explicit per-key bounds, with the same ceiling for union and intersection");
		require(!retainUnmeasured || caps!=null && ceiling>0,
			"Unmeasured-key fallback needs the same explicit family ceiling as its measured neighbors");
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
				if(caps.left[index]<0){
					unavailableSeeds++;
					if(!retainUnmeasured || caps.siblingOnly!=null && caps.siblingOnly[index]){continue;}
					final int fallback=Math.min(padding, ceiling);
					start=(long)center-fallback;
					stop=(long)center+fallback+K;
				}else{
					final long keyStart=(long)center-K/2;
					start=keyStart-caps.left[index];
					stop=keyStart+K-1+caps.right[index];
				}
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
	/** Reconstruct only supplied scanner hits; no second genome scan and no per-base allocation. */
	public static long[] hitKeys(final byte[] bases, final int[] centers, final int k){
		assert(bases!=null && centers!=null && k>0 && k<=31) : "Forward scanner hits retain their sequence and valid fixed seed width";
		final long[] keys=new long[centers.length];
		for(int i=0; i<centers.length; i++){
			final int start=centers[i]-k/2;
			assert(start>=0 && (long)start+k<=bases.length) : "Each supplied seed center owns a complete scanned footprint";
			long word=0;
			for(int j=0; j<k; j++){
				final int value=AminoAcid.baseToNumber[bases[start+j]];
				assert(value>=0) : "Scanner seed hits cannot cross ambiguous sequence";
				word=(word<<2)|value;
			}
			keys[i]=word;
		}
		return keys;
	}
	/** Claim subtraction may shorten geometry but cannot discard its protection metadata. */
	public static int[] part(final int[] window, final int start, final int stop){
		assert(window!=null && window.length>=2 && start>=window[0] && stop<=window[1] && stop>=start)
			: "A surviving claim remainder stays inside its original protected interval";
		final int[] result=window.clone();
		result[0]=start;
		result[1]=stop;
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
			return read(path, member, expectedKeys, capRule(percentile));
		}
		/** Generic tables declare one exact percentile/padding rule and their own positive chain ceiling. */
		public static Caps loadGeneric(final String path, final String member, final int expectedKeys){
			final Caps caps=read(path, member, expectedKeys, null);
			require(caps.chainCeiling>0, "Generic per-key windows require explicit positive chain-ceiling metadata: "+member);
			return caps;
		}
		/** A null required rule is inferred once from the first observed key, never independently per row. */
		private static Caps read(final String path, final String member, final int expectedKeys, final String requiredRule){
			require(path!=null && member!=null && expectedKeys>0, "Per-key bounds need an explicit member, table and fixed-axis size");
			String rule=requiredRule;
			final LongList words=new LongList(expectedKeys);
			final IntList left=new IntList(expectedKeys), right=new IntList(expectedKeys);
			final LineParser1 p=new LineParser1('\t');
			int ceiling=0;
			String ceilingRule=null, population=null, combination=null;
			long line=1;
			try(Input input=new Input(path)){
				RrnaResourceIO.header(input.reader, CAP_HEADER, path);
				for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
					line++;
					p.set(row);
					if(p.termEquals(CEILING_TAG, 0)){
						require(words.size==0 && ceilingRule==null && p.terms()==3,
							"One chain-ceiling metadata row may precede the fixed seed axis");
						ceiling=p.parseInt(1);
						ceilingRule=p.parseString(2);
						validateCeiling(ceiling, ceilingRule);
						continue;
					}
					if(p.termEquals(POPULATION_TAG, 0)){
						require(words.size==0 && population==null && p.terms()==2,
							"One declared construction population may precede the fixed seed axis");
						population=p.parseString(1);
						validatePopulation(population);
						continue;
					}
					if(p.termEquals(COMBINATION_TAG, 0)){
						require(requiredRule==null && words.size==0 && combination==null && p.terms()==2,
							"One generic combining-rule row may precede the seed axis; legacy explicit-policy readers reject it");
						combination=p.parseString(1);
						validateCombination(combination);
						continue;
					}
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
						if(rule==null){
							rule=p.parseString(7);
							validateCapRule(rule);
						}
						require(a>=0 && z>=0 && p.parseLong(6)==K+(long)a+z && p.termEquals(rule, 7),
							"All observed rows require the same declared full-member-capped percentile/padding rule: "+rule);
						left.add(a);
						right.add(z);
					}
					words.add(word);
				}
			}catch(RuntimeException failure){throw new IllegalArgumentException("Invalid per-key bounds at line "+line+" in "+path, failure);}
			require(words.size==expectedKeys, "Every fixed seed key must retain one available or explicitly unavailable bounds row");
			return new Caps(Arrays.copyOf(words.array, words.size), Arrays.copyOf(left.array, left.size),
				Arrays.copyOf(right.array, right.size), ceiling, ceilingRule, population, combination);
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
			this(keys, left_, right_, 0, null);
		}
		/** Legacy tables have no metadata; only their existing external policy may supply a ceiling. */
		public Caps(final long[] keys, final int[] left_, final int[] right_, final int ceiling_, final String rule_){
			this(keys, left_, right_, ceiling_, rule_, null);
		}
		/** Population text is a caller declaration, not independent proof of how genes were extracted. */
		public Caps(final long[] keys, final int[] left_, final int[] right_, final int ceiling_, final String rule_, final String population_){
			this(keys, left_, right_, ceiling_, rule_, population_, null);
		}
		/** An omitted rule preserves union; an explicit generic rule accompanies a positive chain ceiling. */
		public Caps(final long[] keys, final int[] left_, final int[] right_, final int ceiling_, final String rule_,
			final String population_, final String combination_){
			require(keys!=null && left_!=null && right_!=null && keys.length>0
				&& keys.length==left_.length && keys.length==right_.length,
				"Per-key margins require one left/right row on the same nonempty fixed forward seed axis");
			left=left_.clone();
			right=right_.clone();
			require(ceiling_==0 && rule_==null || ceiling_>0 && rule_!=null,
				"Chain-ceiling metadata is either absent as a pair or an explicit positive bound and rule");
			if(rule_!=null){validateCeiling(ceiling_, rule_);}
			chainCeiling=ceiling_;
			chainCeilingRule=rule_;
			if(population_!=null){
				validatePopulation(population_);
				require(ceiling_>0, "Declared construction populations accompany the new explicit chain-ceiling schema");
			}
			constructionPopulation=population_;
			if(combination_!=null){
				validateCombination(combination_);
				require(ceiling_>0, "An explicit generic combining rule requires its protected chain ceiling");
			}
			union=combination_==null || UNION.equals(combination_);
			index=new LongIntMap(keys.length*2);
			int unmeasured=0;
			for(int i=0; i<keys.length; i++){
				require(keys[i]>=0 && keys[i]<(1L<<(2*K)) && index.get(keys[i])<0,
					"The fixed forward seventeen-mer axis must contain unique valid words");
				require(left[i]>=-1 && right[i]>=-1 && (left[i]<0)==(right[i]<0),
					"Unavailable margins use -1 at both ends; observed margins must be nonnegative");
				index.put(keys[i], i);
				if(left[i]<0){unmeasured++;}
			}
			unmeasuredKeys=unmeasured;
			siblingOnly=null;
			siblingOnlyKeys=0;
		}
		/** Shares immutable margins and indexing; only the new attribution mask is owned here. */
		private Caps(final Caps source, final boolean[] excluded, final int count){
			assert(excluded!=null && excluded.length==source.left.length && count>0 && count<=source.unmeasuredKeys)
				: "Sibling exclusion can mark only existing unobserved member rows before worker publication";
			index=source.index;
			left=source.left;
			right=source.right;
			chainCeiling=source.chainCeiling;
			chainCeilingRule=source.chainCeilingRule;
			constructionPopulation=source.constructionPopulation;
			union=source.union;
			unmeasuredKeys=source.unmeasuredKeys;
			siblingOnly=excluded;
			siblingOnlyKeys=count;
		}
		/** Cold group binding excludes a key only with positive observations in a sibling.
		 * Unseen-in-all keys remain eligible for the ordinary capped fallback. */
		Caps withoutSiblingOnly(final Caps sibling){
			require(sibling!=null, "Member attribution needs an actual loaded sibling table, not inferred absence");
			final long[] keys=index.keys();
			final int[] rows=index.values();
			assert(keys.length==rows.length) : "LongIntMap exposes parallel immutable key/value arrays for the bounds row index";
			boolean[] excluded=null;
			int count=siblingOnlyKeys;
			for(int i=0; i<keys.length; i++){
				if(keys[i]<0){continue;}
				final int row=rows[i];
				if(left[row]>=0 || sibling.leftMargin(keys[i])<0 || siblingOnly!=null && siblingOnly[row]){continue;}
				if(excluded==null){excluded=siblingOnly==null ? new boolean[left.length] : siblingOnly.clone();}
				excluded[row]=true;
				count++;
			}
			return excluded==null ? this : new Caps(this, excluded, count);
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
		private final boolean[] siblingOnly;
		public final int chainCeiling;
		public final String chainCeilingRule;
		public final String constructionPopulation;
		/** Combines distinct-key envelopes only; repeated-key envelopes and across-member pooling remain unions. */
		public final boolean union;
		public final int unmeasuredKeys;
		/** Unobserved here but observed in a loaded sibling; these keys cannot open a window. */
		public final int siblingOnlyKeys;
	}
	private static void validateCeiling(final int ceiling, final String rule){
		require(ceiling>0 && (CEILING_P999.equals(rule) || CEILING_MAX.equals(rule)),
			"Chain ceilings require a positive gene-length-minus17 value and the declared sample-size rule");
	}
	private static void validatePopulation(final String population){
		require(POPULATION_ORIGINALS.equals(population) || POPULATION_SCORED.equals(population),
			"Bounds must explicitly declare original genes or scored reference intervals in-corpus");
	}
	private static void validateCombination(final String rule){
		require(UNION.equals(rule) || INTERSECTION.equals(rule), "Generic bounds combine distinct keys with exactly UNION or INTERSECTION");
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
	public static final String CEILING_TAG="#window_ceiling_bp";
	public static final String CEILING_P999="P999_MINUS17_N_GE_1000";
	public static final String CEILING_MAX="MAX_MINUS17_N_LT_1000";
	public static final String POPULATION_TAG="#construction_population";
	public static final String POPULATION_ORIGINALS="original genes";
	public static final String POPULATION_SCORED="built from scored reference intervals, in-corpus";
	public static final String COMBINATION_TAG="#window_combine";
	public static final String UNION="UNION", INTERSECTION="INTERSECTION";
	/** Legacy explicit readers retain their integer-percentile plus-ten contract. */
	public static String capRule(final int percentile){
		return capRule((double)percentile, 10);
	}
	/** Shared producer/consumer label records the actual nearest rank and nonnegative end padding. */
	public static String capRule(final double percentile, final int padding){
		require((percentile==90 || percentile==95 || percentile==98 || percentile==99 || percentile==99.9) && padding>=0,
			"Bounds require percentile90/95/98/99/99.9 and nonnegative added bases");
		if(percentile==99 && padding==10){return CAP_RULE;}
		final String text=percentile==99.9 ? "99.9" : Integer.toString((int)percentile);
		return "P"+text+"_PLUS"+padding+"_MEMBER_MAX_MINUS17";
	}
	/** Parse one cold metadata value, rejecting legacy diagnostics, noncanonical text and unknown rules. */
	private static void validateCapRule(final String rule){
		final String suffix="_MEMBER_MAX_MINUS17";
		final int plus=rule.indexOf("_PLUS"), end=rule.length()-suffix.length();
		require(rule.startsWith("P") && rule.endsWith(suffix) && plus>1 && plus+5<end,
			"Generic bounds require an explicit full-member-capped percentile and added-base estimator: "+rule);
		final double percentile=Double.parseDouble(rule.substring(1, plus));
		final int padding=Integer.parseInt(rule.substring(plus+5, end));
		require(rule.equals(capRule(percentile, padding)), "Generic bounds estimator text must be canonical: "+rule);
	}
	private static final Comparator<int[]> BY_SPAN=(a, b)->{
		final int start=Integer.compare(a[0], b[0]);
		return start==0 ? Integer.compare(a[1], b[1]) : start;
	};
}
