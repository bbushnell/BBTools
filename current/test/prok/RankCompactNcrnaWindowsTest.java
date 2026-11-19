package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import cardinality.DynamicDemiLog;
import idaligner.Tracer;
import map.IntHashMap;
import map.IntHashSet;
import parse.LineParser1;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.BuildNcrnaDdlProfiles;
import prok.VerifyNcrnaDdlProfiles;
import prok.RankCompactNcrnaWindows;
import prok.RankDualNcrnaWindows;
import prok.NcrnaDdlQueryKeys;
import prok.NcrnaDdlKeyFormat;
import prok.NcrnaDdlLayoutKeys;
import prok.NcrnaDdlLayoutVotes;
import prok.NcrnaDdlLayoutModels;
import prok.NcrnaDdlMembers;
import prok.NcrnaDdlModelProfiles;
import prok.NcrnaDdlPartitionKeys;
import prok.NcrnaDdlProfile;
import prok.NcrnaDdlSavedData;
import shared.Tools;
import shared.Shared;
import structures.ByteBuilder;

/** Uses the original two-model saved-cache witness, never a genome or new oracle.
 * Literal report counts, known ranks, key sizes and original UTF-8 rows are
 * checked with alignment forbidden. The separate vote fixture tests collisions.
 * @author Brian, Raiden
 */
public final class RankCompactNcrnaWindowsTest{

	static void run(final RankSavedNcrnaWindowsTest.Fixture f){
		consensusCollisions();
		keyFormats();
		layoutKeys();
		layoutVotes();
		final RankDualNcrnaWindows.Config c=RankDualNcrnaWindowsTest.config(f, "compact", 1);
		c.identities=f.dir+"/one";
		final boolean global=Tracer.GLOBAL;
		Tracer.GLOBAL=true;
		try{
			RankCompactNcrnaWindows.run(c);
			DynamicDemiLog.setExponent(6);
			RankCompactNcrnaWindows.main(new String[]{"panels="+c.panels, "selection="+c.selection,
				"baseline="+c.baseline, "identities="+c.identities, "members="+c.members,
				"ledger="+c.ledger, "assignments="+c.assignments, "census="+c.census, "summary="+c.summary,
				"ref="+c.references, "out="+c.output+"_cli", "expectedmodels=2", "expectedgroups=3", "t=1"});
			check(DynamicDemiLog.exponentBits()==6, "Compact CLI must restore its embedding caller's native exponent");
		}finally{
			Tracer.GLOBAL=global;
			DynamicDemiLog.setExponent(4);
		}
		for(final String suffix:SUFFIXES){equal(c.output+suffix, c.output+"_cli"+suffix);}
		final LineParser1 p=new LineParser1('\t');
		int count=0;
		try(Input in=new Input(c.output+".methods.tsv")){
			check(in.reader.nextLine()!=null, "Compact method definitions must precede all rows");
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				check(p.terms()==15 && p.parseLong(2)==2L*p.parseInt(1) && p.parseLong(3)==4L*p.parseInt(1)
					&& p.parseInt(4)==2, "Key payload uses exactly two bytes per bin per original model, excluding native implicit keys/postings");
				if(p.termStartsWith("consensus_weight1000000_", 0)){
					check(p.parseDouble(5)>=1 && p.parseInt(6)>=1 && p.parseInt(7)>=p.parseInt(6)
						&& p.parseDouble(8)>0 && p.parseInt(9)>0 && p.parseInt(10)>=p.parseInt(9),
						"Always-win rows require consensus occupancy and length mean/range, not only a rank");
				}
				count++;
			}
		}
		check(count==22, "Keep two accepted controls, eight methods each at512/1024 and four at2048; no4096study");
		count=0;
		try(Input in=new Input(c.output+".profile_models.tsv")){
			check(in.reader.nextLine()!=null, "Per-model occupancy/length evidence needs a header");
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				check(p.terms()==13 && p.parseInt(1)>=0 && p.parseInt(1)<2 && p.parseInt(3)>0
					&& p.parseInt(4)==1 && p.parseInt(6)<=p.parseInt(5), "Every compact method retains real membership and bounded consensus matches");
				check(p.parseInt(9)>=p.parseInt(5) && p.parseInt(10)<=p.parseInt(5)
					&& p.parseInt(11)==0 && p.termEquals("NA", 12),
					"Raw-word collision counts precede key compression; million and historical always keys agree in the literal fixture");
				if(p.termStartsWith("consensus_weight1000000_", 0)){
					check(p.parseInt(5)==p.parseInt(6), "Consensus must own every occupied consensus bin under always-win construction");
				}
				count++;
			}
		}
		check(count==40, "Every one of twenty compact methods must retain both original model rows");
		count=0;
		try(Input in=new Input(c.output+".summary.tsv")){
			check(in.reader.nextLine()!=null, "Summary needs explicit original membership and identity bins");
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				check(p.terms()==15, "Every summary row retains both last-input-loss controls and storage payload");
				if(p.termEquals("NOT_SEEN", 1) && (p.termEquals("ALL", 2) || p.termEquals("Q1_LOW", 2))){
					check(p.parseInt(4)==1 && p.parseInt(5)==1 && p.parseInt(6)==1 && p.parseInt(7)==1
						&& p.parseInt(11)==0 && p.parseInt(13)==0,
						"The literal A-gene common locus ranks first in every method and preserves its only reachable input");
				}
				count++;
			}
		}
		check(count==396, "All twenty-two methods keep three original membership states and six identity-bin views");
		check(rows(c.output+".inputs.tsv")==4 && rows(c.output+".loci.tsv.gz")==133
			&& rows(c.output+".last_losses.tsv")==1,
			"Zero-target inputs and all unrankable/unknown parents survive; no fixture input loses its final locus");
		boolean rejected=false;
		try{RankCompactNcrnaWindows.run(c);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Compact reports cannot overwrite accepted output");
		c.output=f.dir+"/compact_missing_identity";
		c.identities=null;
		rejected=false;
		try{RankCompactNcrnaWindows.run(c);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected && !new File(c.output+".methods.tsv").exists(), "Missing saved identities must fail instead of measuring new ones");
		partitions(f);
		fineWeights(f);
		System.err.println("NCRNA_COMPACT_S1_SAVED_RANKS_TEST_PASS");
	}

	private static void partitions(final RankSavedNcrnaWindowsTest.Fixture f){
		final RankDualNcrnaWindows.Config c=RankDualNcrnaWindowsTest.config(f, "compact_partitions", 1);
		c.identities=f.dir+"/one";
		partitionKeys(c);
		final boolean global=Tracer.GLOBAL;
		Tracer.GLOBAL=true;
		try{
			RankCompactNcrnaWindows.runPartitions(c);
			RankCompactNcrnaWindows.main(new String[]{"panels="+c.panels, "selection="+c.selection,
				"baseline="+c.baseline, "identities="+c.identities, "members="+c.members,
				"ledger="+c.ledger, "assignments="+c.assignments, "census="+c.census, "summary="+c.summary,
				"ref="+c.references, "out="+c.output+"_cli", "expectedmodels=2", "expectedgroups=3", "partitions=t", "t=1"});
		}finally{Tracer.GLOBAL=global;}
		for(final String suffix:SUFFIXES){equal(c.output+suffix, c.output+"_cli"+suffix);}
		check(rows(c.output+".methods.tsv")==12 && rows(c.output+".profile_models.tsv")==19,
			"Only nine partition geometries and two saved controls run; accepted S1 profiles are not rebuilt");
		check(rows(c.output+".summary.tsv")==199 && rows(c.output+".loci.tsv.gz")==67
			&& rows(c.output+".inputs.tsv")==4, "All membership/identity views, unrankable parents and zero-target inputs remain");
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(c.output+".profile_models.tsv")){
			check(input.reader.nextLine()!=null, "Partition model statistics require a header");
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				final long bytes=p.parseLong(8);
				check((bytes==1024 || bytes==2048 || bytes==4096) && p.parseInt(4)==1
					&& p.parseInt(5)==p.parseInt(6) && p.parseInt(5)<=bytes/8
					&& p.parseInt(7)>=p.parseInt(5),
					"Consensus owns only its reserved quarter, with exact uint16 width and original member count");
			}
		}
		try(Input input=new Input(c.output+".summary.tsv")){
			check(input.reader.nextLine()!=null, "Partition ranks retain the same accepted identity strata");
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				if(p.termEquals("NOT_SEEN", 1) && (p.termEquals("ALL", 2) || p.termEquals("Q1_LOW", 2))){
					check(p.parseInt(4)==1 && p.parseInt(5)==1 && p.parseInt(6)==1 && p.parseInt(7)==1,
						"The literal A-gene witness ranks first without any new identity or alignment work");
				}
			}
		}
		resourceParity(c);
		System.err.println("NCRNA_COMPACT_PARTITIONS_SAVED_RANKS_TEST_PASS");
	}

	/** Literal span enumeration checks plain construction keys and exact queries.
	 * The new consensus-only option must leave all member areas unchanged. */
	private static void partitionKeys(final RankDualNcrnaWindows.Config c){
		final byte[] sequence=Prok16sPartialRescueTest.random(900, 71523);
		sequence[180]='N';
		sequence[510]='R';
		for(final int width:new int[]{512,1024}){
			for(final boolean both:new boolean[]{true,false}){
				final NcrnaDdlPartitionKeys codec=new NcrnaDdlPartitionKeys(7, width, both);
				final char[] expected=new char[width];
				final IntHashSet queries=new IntHashSet(4096);
				literalPartition(sequence, width, both, expected, queries);
				check(Arrays.equals(expected, codec.consensusKeys(sequence)),
					"Plain consensus maxima must equal independent complete-span enumeration");
				codec.build(sequence);
				check(codec.keyCount()==queries.size(), "Consensus-only queries cannot retain hidden short probes or lose member probes");
				for(int i=0; i<codec.keyCount(); i++){
					check(queries.contains(codec.packedKey(i)), "Every raw query key must belong to its exact declared member/consensus area");
				}
				final char[] queryMaxima=new char[width];
				for(final int packed:queries.toArray()){
					final int bin=packed>>>16;
					if((char)packed>queryMaxima[bin]){queryMaxima[bin]=(char)packed;}
				}
				codec.build(sequence, true);
				int next=0;
				for(int bin=0; bin<width; bin++){
					if(queryMaxima[bin]==0){continue;}
					check(next<codec.keyCount() && codec.packedKey(next++)==((bin<<16)|queryMaxima[bin]),
						"Legacy query-DDL probes must retain the independent raw-set maximum of each nonempty bin");
				}
				check(next==codec.keyCount(), "Legacy query-DDL cannot retain displaced raw keys or empty-bin matches");
				for(final char key:codec.consensusKeys(new byte[]{'N','R'})){
					check(key==0, "Undefined-only replacement must clear every old consensus key");
				}
			}
		}
		final ArrayList<String> names=new ArrayList<String>();
		final ArrayList<byte[]> sequences=new ArrayList<byte[]>();
		NcrnaDdlSavedData.sequences(c.references, read->{
			names.add(read.id);
			sequences.add(read.bases);
		});
		final byte[][] consensus=sequences.toArray(new byte[sequences.size()][]);
		final NcrnaDdlMembers members=NcrnaDdlMembers.load(c.members, c.ledger, c.assignments,
			c.census, c.summary, names.toArray(new String[names.size()]), consensus);
		memberLayoutGrid(members);
		final NcrnaDdlProfile[] legacy=NcrnaDdlModelProfiles.buildPartition(members, consensus, 7, 512);
		final NcrnaDdlProfile[] both=NcrnaDdlModelProfiles.buildPartition(members, consensus, 7, 512, true);
		final NcrnaDdlProfile[] only=NcrnaDdlModelProfiles.buildPartition(members, consensus, 7, 512, false);
		for(int model=0; model<members.models(); model++){
			check(Arrays.equals(legacy[model].keys(), both[model].keys()) && both[model].members==legacy[model].members,
				"Explicit default consensus content must retain accepted keys and original member counts");
			final NcrnaDdlPartitionKeys longCodec=new NcrnaDdlPartitionKeys(7, 512, false);
			final char[] longKeys=longCodec.consensusKeys(consensus[model]);
			for(int bucket=0; bucket<512; bucket++){
				if(bucket<384){check(both[model].keys()[bucket]==only[model].keys()[bucket], "Consensus options cannot change a member partition");}
				else{
					check(only[model].keys()[bucket]==longKeys[bucket],
						"Long-only consensus bins must retain exactly the independent masked-ten-base maxima");
				}
			}
			for(int previous=0; previous<model; previous++){
				check(both[model].keys()!=both[previous].keys() && only[model].keys()!=only[previous].keys(),
					"Each model must own its plain keys after construction scratch is reused");
			}
		}
		System.err.println("NCRNA_DDL_PLAIN_PARTITION_KEYS_TEST_PASS");
	}

	/** Check every leading-zero tier against an independent bit reader, preserving
	 * native exponent masking and Java's masked shift distance. */
	private static void keyFormats(){
		for(final int[] fields:new int[][]{{6,10},{5,10},{4,11},{4,12}}){
			final NcrnaDdlKeyFormat format=new NcrnaDdlKeyFormat(fields[0], fields[1]);
			check(format.encode(0)==(1<<fields[1])-1 && format.encode(-1L)==1,
				"Native exponent64 wraps to zero while an all-one hash reserves key zero for empty");
			for(int zeros=0; zeros<64; zeros++){
				final long leading=1L<<(63-zeros);
				for(final long suffix:new long[]{0L, -1L, 0x5555555555555555L, 0x2aaaaaaaaaaaaaaaL}){
					final long hash=leading|(suffix&(leading-1));
					final int expected=literalKey(hash, fields[0], fields[1]);
					check(format.encode(hash)==expected, "Explicit exponent/mantissa fields must match independently read bits at every exponent");
				}
			}
		}
		System.err.println("NCRNA_DDL_EXPLICIT_KEY_FORMATS_TEST_PASS");
	}

	private static int literalKey(final long hash, final int exponent, final int mantissa){
		assert(exponent>=1 && exponent<=6 && mantissa>=1 && exponent+mantissa<=16)
			: "The independent bit reader must use the same representable absolute-key fields as its declared fixture";
		int zeros=0;
		while(zeros<64 && ((hash>>>(63-zeros))&1)==0){zeros++;}
		zeros%=1<<exponent;
		final int shift=(63-mantissa-zeros)&63;
		int fraction=0;
		for(int i=0; i<mantissa; i++){
			final int position=shift+i;
			final int bit=position>=64 ? 0 : (int)((hash>>>position)&1);
			fraction|=(1-bit)<<i;
		}
		return Math.max(1, zeros*(1<<mantissa)+fraction);
	}

	/** Separate span enumeration checks full-width plain keys, explicit lengths,
	 * masks, partition splits, query forms and reused empty scratch. */
	private static void layoutKeys(){
		final byte[] sequence="ACGTAGCAGTCAGTACGTTGNNacgtagctagctagcatgcRACGTAGCAGTACGTACGT".getBytes(StandardCharsets.US_ASCII);
		for(final int shortK:new int[]{6,7,8}){
			for(final int longK:new int[]{9,10,11,12}){
				for(final boolean masked:new boolean[]{false,true}){
					if(!masked && longK!=10){continue;}
					final NcrnaDdlLayoutKeys.Raw raw=new NcrnaDdlLayoutKeys.Raw(shortK, longK, masked);
					for(final int bins:new int[]{128,192,256,384,512}){
						for(final int split:new int[]{0,bins/4,bins/3,bins/2}){
							layoutKeyCase(sequence, raw, bins, split, true);
							if(split==0){layoutKeyCase(sequence, raw, bins, split, false);}
						}
					}
				}
			}
		}
		final NcrnaDdlLayoutKeys.Raw raw=new NcrnaDdlLayoutKeys.Raw();
		raw.build("AAAAAAANCCCCCCCCCC".getBytes(StandardCharsets.US_ASCII), true);
		final long shortC=5461L, maskedLongC=349525L&~(3L<<8);
		check(raw.maxExponent(true)==Math.max(Long.numberOfLeadingZeros(Tools.hash64shift((1L<<62)^12345L)),
			Long.numberOfLeadingZeros(Tools.hash64shift(shortC^(1L<<62)^12345L)))
			&& raw.maxExponent(false)==Long.numberOfLeadingZeros(Tools.hash64shift(maskedLongC^(1L<<63)^12345L))
			&& raw.legacyMaxExponent(true)==Math.max(Long.numberOfLeadingZeros(Tools.hash64shift(12345L)),
			Long.numberOfLeadingZeros(Tools.hash64shift(shortC^12345L)))
			&& raw.legacyMaxExponent(false)==Long.numberOfLeadingZeros(Tools.hash64shift(maskedLongC^12345L)),
			"Literal A7/C7/C10 words must report actual unwrapped exponents in both declared hash namespaces");
		raw.build(sequence);
		check(raw.legacyMaxExponent(true)==-1 && raw.legacyMaxExponent(false)==-1,
			"Ordinary consensus construction must clear and skip the optional legacy diagnostic hash");
		System.err.println("NCRNA_DDL_LAYOUT_KEYS_TEST_PASS");
	}

	private static void layoutKeyCase(final byte[] sequence, final NcrnaDdlLayoutKeys.Raw raw,
			final int bins, final int split, final boolean both){
		assert(both || split==0) : "Long-only consensus fixtures use one undivided key array";
		if(Integer.bitCount(split==0 ? bins : split)!=1 || split>0 && Integer.bitCount(bins-split)!=1){
			boolean rejected=false;
			try{new NcrnaDdlLayoutKeys(bins, split, both, new NcrnaDdlKeyFormat(4, 12), raw.shortK, raw.longK, raw.maskedLong);}
			catch(IllegalArgumentException expected){rejected=true;}
			check(rejected, "A native layout must reject a non-power-of-two actual part");
			return;
		}
		final NcrnaDdlLayoutKeys query=new NcrnaDdlLayoutKeys(bins, split, both,
			new NcrnaDdlKeyFormat(4, 12), raw.shortK, raw.longK, raw.maskedLong);
		final IntHashSet expected=new IntHashSet(256);
		final char[] maxima=new char[bins];
		literalLayout(sequence, bins, split, both, raw.shortK, raw.longK, raw.maskedLong, expected, maxima);
		raw.build(sequence);
		query.build(raw);
		check(Arrays.equals(maxima, query.maximumKeys()), "Plain maxima must retain exactly the greatest compressed key in each bin");
		int occupied=0;
		for(final char key:maxima){if(key!=0){occupied++;}}
		check(query.keyCount()==occupied, "The native voting prefix retains one greatest nonempty key per bin");
		final boolean[] seenBins=new boolean[bins];
		for(int i=0; i<query.keyCount(); i++){
			final int packed=query.packedKey(i), bin=packed>>>16;
			check(!seenBins[bin] && maxima[bin]!=0 && maxima[bin]==(char)packed,
				"A query DDL cannot emit zero, repeat a bin or retain a displaced raw key");
			seenBins[bin]=true;
		}
		raw.build(new byte[]{'N','R'});
		query.build(raw);
		check(query.keyCount()==0
			&& raw.maxExponent(true)==-1 && raw.maxExponent(false)==-1,
			"Undefined-only replacement must clear raw/query-DDL prefixes and observed exponents");
		for(final char key:query.maximumKeys()){check(key==0, "Empty replacement clears old layout maxima");}
	}

	private static void literalLayout(final byte[] sequence, final int bins, final int shortBins,
			final boolean both, final int shortK, final int longK, final boolean masked,
			final IntHashSet queries, final char[] maxima){
		assert(sequence!=null && shortK>0 && shortK<longK && longK<=15 && shortBins>=0 && shortBins<bins
			&& queries!=null && maxima!=null && maxima.length==bins)
			: "Literal span enumeration needs complete fixture bytes and the declared key-only projection geometry";
		for(final int k:new int[]{shortK,longK}){
			if(k==shortK && !both){continue;}
			for(int start=0; start+k<=sequence.length; start++){
				long word=0;
				boolean valid=true;
				for(int p=start; p<start+k; p++){
					final int code="ACTG".indexOf(Character.toUpperCase((char)sequence[p]));
					if(code<0){
						valid=false;
						break;
					}
					word=4*word+code;
				}
				if(!valid){continue;}
				if(k==longK && masked){word&=~(3L<<(2*((longK-1)/2)));}
				final long tag=shortBins>0 ? 0 : k==shortK ? 1L<<62 : 1L<<63;
				final long hash=Tools.hash64shift(word^tag^12345L);
				final int count=shortBins>0 ? (k==shortK ? shortBins : bins-shortBins) : bins;
				final int bin=(int)(hash&(count-1))+(shortBins>0 && k==longK ? shortBins : 0);
				final char key=(char)literalKey(hash, 4, 12);
				queries.add((bin<<16)|key);
				if(key>maxima[bin]){maxima[bin]=key;}
			}
		}
	}

	/** Compressed-key equality defines a vote candidate, regardless of source
	 * length. Bonuses do not change real membership or mutate reusable counts. */
	private static void layoutVotes(){
		final NcrnaDdlLayoutVotes ties=new NcrnaDdlLayoutVotes(1);
		ties.add(new char[]{10});
		ties.add(new char[]{11});
		final NcrnaDdlLayoutVotes.Result highTie=ties.finish();
		check(highTie.keys[0]==11 && highTie.members==2, "Greatest full key must break equal-support plain-key ties");
		ties.add(new char[]{10});
		final NcrnaDdlLayoutVotes.Result plurality=ties.finish();
		check(plurality.keys[0]==10 && plurality.members==3 && highTie.keys[0]==11 && highTie.members==2,
			"Equal compressed keys combine support, outvote a greater singleton and leave earlier results independent");
		final NcrnaDdlLayoutKeys.Raw raw=new NcrnaDdlLayoutKeys.Raw();
		raw.build("AAAAAAA".getBytes(StandardCharsets.US_ASCII));
		final NcrnaDdlLayoutKeys consensus=new NcrnaDdlLayoutKeys(1, false, true, new NcrnaDdlKeyFormat(4, 12));
		consensus.buildWords(new int[]{0,0}, 2, new int[0], 0);
		check(consensus.keyCount()==1, "Repeated consensus words receive the bonus only once per native winning bin/key");
		final char key=consensus.maximumKeys()[0];
		final char other=(char)(key==65535 ? key-1 : key+1);
		final NcrnaDdlLayoutVotes votes=new NcrnaDdlLayoutVotes(1);
		votes.add(new char[]{other});
		final NcrnaDdlLayoutVotes.Result noBonus=votes.finish(consensus, 0);
		final NcrnaDdlLayoutVotes.Result bonus=votes.finish(consensus, 2);
		check(noBonus.keys[0]==other && noBonus.members==1 && bonus.keys[0]==key && bonus.members==1
			&& votes.finish().keys[0]==other && noBonus.keys[0]==other,
			"A consensus-only candidate can win the pooled bonus without inventing real voters or changing prior results");
		votes.add(new char[]{key});
		check(votes.finish(consensus, Integer.MAX_VALUE).keys[0]==key,
			"Bonus arithmetic must remain representable when real support is added to the largest declared integer bonus");
		System.err.println("NCRNA_DDL_LAYOUT_VOTES_TEST_PASS");
	}

	/** The two-model, one-member-per-model witness exposes scratch aliasing:
	 * each geometry's final model must equal that member's independently built
	 * maximum sketch after every other model/geometry has reused the scratch. */
	private static void memberLayoutGrid(final NcrnaDdlMembers members){
		final NcrnaDdlLayoutKeys[] layouts=new NcrnaDdlLayoutKeys[6];
		final int[] bins={256,192,256,384,512,512}, splits={0,64,128,128,0,256};
		for(int i=0; i<layouts.length; i++){
			layouts[i]=new NcrnaDdlLayoutKeys(bins[i], splits[i], true,
				new NcrnaDdlKeyFormat(4, 12), 7, 10, true);
		}
		final NcrnaDdlLayoutModels.Result grid=NcrnaDdlLayoutModels.buildMembers(members, layouts);
		check(grid.voters==members.keptMembers() && grid.voters==2,
			"The complete grid retains exactly the toy's two original voters, not one invented population per geometry");
		final NcrnaDdlLayoutKeys.Raw raw=new NcrnaDdlLayoutKeys.Raw();
		for(int geometry=0; geometry<layouts.length; geometry++){
			for(int model=0; model<members.models(); model++){
				check(members.members(model).length==1, "The literal grid witness needs one original member per model");
				raw.build(members.members(model)[0]);
				layouts[geometry].build(raw);
				final NcrnaDdlLayoutVotes.Result profile=grid.profiles[geometry][model];
				check(profile.members==1 && Arrays.equals(profile.keys, layouts[geometry].maximumKeys()),
					"Each retained model/geometry must own its exact selected plain keys after all scratch reuse");
				for(int other=0; other<geometry; other++){
					check(profile.keys!=grid.profiles[other][model].keys,
						"A geometry result must not borrow another geometry's mutable profile storage");
				}
			}
		}
		System.err.println("NCRNA_DDL_LAYOUT_MODEL_GRID_TEST_PASS");
	}

	/** New weights only: accepted 10/100 rows remain external scientific controls. */
	private static void fineWeights(final RankSavedNcrnaWindowsTest.Fixture f){
		final RankDualNcrnaWindows.Config c=RankDualNcrnaWindowsTest.config(f, "compact_fine", 1);
		c.identities=f.dir+"/one";
		final boolean global=Tracer.GLOBAL;
		Tracer.GLOBAL=true;
		try{
			RankCompactNcrnaWindows.runFineWeights(c);
			RankCompactNcrnaWindows.main(new String[]{"panels="+c.panels, "selection="+c.selection,
				"baseline="+c.baseline, "identities="+c.identities, "members="+c.members,
				"ledger="+c.ledger, "assignments="+c.assignments, "census="+c.census, "summary="+c.summary,
				"ref="+c.references, "out="+c.output+"_cli", "expectedmodels=2", "expectedgroups=3", "fineweights=t", "t=1"});
		}finally{Tracer.GLOBAL=global;}
		for(final String suffix:SUFFIXES){equal(c.output+suffix, c.output+"_cli"+suffix);}
		check(rows(c.output+".methods.tsv")==9 && rows(c.output+".profile_models.tsv")==13
			&& rows(c.output+".loci.tsv.gz")==49 && rows(c.output+".summary.tsv")==145,
			"Only six new weights plus two saved controls retain every original parent and reporting stratum");
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(c.output+".methods.tsv")){
			check(input.reader.nextLine()!=null && input.reader.nextLine()!=null && input.reader.nextLine()!=null,
				"The fine report retains its header and two accepted saved control rows");
			for(final int weight:new int[]{2,3,5,8,20,50}){
				p.set(input.reader.nextLine());
				check(p.termEquals("consensus_weight"+weight+"_1024", 0) && p.parseInt(1)==1024,
					"The fine pass must not silently rerun accepted weights or introduce another width");
			}
			check(input.reader.nextLine()==null, "No unrequested fine-weight method may be appended");
		}
		System.err.println("NCRNA_COMPACT_FINE_WEIGHTS_TEST_PASS");
	}

	/** An explicit k-by-position loop is independent of the production rolling loop. */
	private static void literalPartition(final byte[] sequence, final int width, final boolean both,
			final char[] winners, final IntHashSet queries){
		final int quarter=width/4;
		for(final int k:new int[]{7,10}){
			for(int start=0; start+k<=sequence.length; start++){
				long word=0;
				boolean valid=true;
				for(int i=start; i<start+k; i++){
					final byte base=sequence[i];
					if(base!='A' && base!='C' && base!='G' && base!='T'){
						valid=false;
						break;
					}
					word=(word<<2)|((base>>1)&3);
				}
				if(!valid){continue;}
				if(k==10){word&=~(3L<<8);}
				final long member=Tools.hash64shift(word^12345L);
				final int memberBucket=k==7 ? (int)(member&(quarter-1)) : quarter+(int)(member&(2*quarter-1));
				queries.add((memberBucket<<16)|literalKey(member, 4, 12));
				if(k==7 && !both){continue;}
				final long hash=Tools.hash64shift(word^(k==7 ? 1L<<62 : 1L<<63)^12345L);
				final int offset=(int)(hash&(quarter-1)), bucket=3*quarter+offset;
				final char key=(char)literalKey(hash, 4, 12);
				queries.add((bucket<<16)|key);
				if(key>winners[bucket]){winners[bucket]=key;}
			}
		}
	}

	/** The actual builder/loader/scorer must reproduce the existing saved report,
	 * including UTF-8 panel names and every unavailable parent. No aligner runs. */
	private static void resourceParity(final RankDualNcrnaWindows.Config saved){
		final BuildNcrnaDdlProfiles.Config build=new BuildNcrnaDdlProfiles.Config();
		build.fasta=saved.members;
		build.references=saved.references;
		build.ledger=saved.ledger;
		build.assignments=saved.assignments;
		build.census=saved.census;
		build.summary=saved.summary;
		build.models=2;
		build.members=2;
		final String left=saved.output+"_profiles512.tsv", right=saved.output+"_profiles1024.tsv";
		final boolean global=Tracer.GLOBAL;
		Tracer.GLOBAL=true;
		try{
			build.partitionBins=512;
			build.output=left;
			BuildNcrnaDdlProfiles.build(build);
			build.partitionBins=1024;
			build.output=right;
			BuildNcrnaDdlProfiles.build(build);
			final VerifyNcrnaDdlProfiles.Config c=parityConfig(saved, left, right, "_parity.tsv");
			VerifyNcrnaDdlProfiles.run(c);
			final LineParser1 p=new LineParser1('\t');
			try(Input input=new Input(c.output)){
				check(input.reader.nextLine()!=null, "A completed parity table requires its count-unit header");
				for(final int width:new int[]{512,1024}){
					p.set(input.reader.nextLine());
					check(p.terms()==8 && p.parseInt(0)==width && p.parseInt(1)==2 && p.parseInt(2)==3
						&& p.parseInt(3)>0 && p.parseInt(4)==6 && p.parseInt(5)==6
						&& p.parseInt(6)+p.parseInt(7)==6 && p.parseInt(7)>0,
						"Both resources must retain all six parents, unavailable ranks and the third zero-target input");
				}
				check(input.reader.nextLine()==null, "Only the two selected caller widths belong in resource parity");
			}
			DynamicDemiLog.setExponent(6);
			final int threads=Shared.threads();
			parityCli(c);
			check(DynamicDemiLog.exponentBits()==6 && Shared.threads()==threads,
				"Ordinary verification and removed-flag rejection must preserve the embedding caller's exponent and thread count");
			final VerifyNcrnaDdlProfiles.Config swapped=parityConfig(saved, right, left, "_parity_swapped.tsv");
			rejectParity(swapped, "Resource slots");
			parityFault(saved, left, right, "rank", "rank differs");
			parityFault(saved, left, right, "input", "panel/input/parent");
			parityFault(saved, left, right, "missing", "omits original parents");
			parityFault(saved, left, right, "duplicate", "extra or duplicate");
		}finally{
			Tracer.GLOBAL=global;
			DynamicDemiLog.setExponent(4);
		}
		System.err.println("NCRNA_DDL_RESOURCE_PARITY_FIXTURE_PASS");
	}

	/** Preserve the ordinary CLI result and reject the retired diagnostic clearly. */
	private static void parityCli(final VerifyNcrnaDdlProfiles.Config c){
		assert(c!=null) : "The ordinary CLI witness retains the complete literal-checked saved fixture";
		final String output=c.output+"_cli.tsv", log=output+".log";
		final ArrayList<String> args=new ArrayList<String>();
		args.add("in="+c.profile512);
		args.add("in2="+c.profile1024);
		args.add("panels="+c.panels);
		args.add("selection="+c.selection);
		args.add("ref="+c.references);
		args.add("ranks="+c.ranks);
		args.add("out="+output);
		args.add("outstream="+log);
		args.add("expectedmodels=2");
		args.add("expectedgroups=3");
		args.add("t=1");
		VerifyNcrnaDdlProfiles.main(args.toArray(new String[args.size()]));
		equal(c.output, output);
		final byte[] normal="NCRNA_DDL_RESOURCE_SAVED_RANK_PARITY_PASS".getBytes(StandardCharsets.US_ASCII);
		final byte[] index="NCRNA_DDL_SMALL_INDEX_SAVED_PARITY_PASS".getBytes(StandardCharsets.US_ASCII);
		int normalMarkers=0, indexMarkers=0;
		try(Input in=new Input(log)){
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				if(Arrays.equals(row, normal)){normalMarkers++;}
				if(Arrays.equals(row, index)){indexMarkers++;}
			}
		}
		check(normalMarkers==1 && indexMarkers==0,
			"Ordinary saved-rank verification must not advertise a comparison of the shared index with itself");
		args.set(args.indexOf("out="+output), "out="+output+".retired");
		args.set(args.indexOf("outstream="+log), "outstream="+log+".retired");
		args.add("indexparity=t");
		boolean rejected=false;
		try{VerifyNcrnaDdlProfiles.main(args.toArray(new String[args.size()]));}
		catch(IllegalArgumentException expected){
			check(expected.getMessage().contains("has been retired"), "The removed option must identify its deliberate retirement");
			rejected=true;
		}
		check(rejected, "Removed diagnostic use must fail rather than silently compare the shared index with itself");
	}

	private static VerifyNcrnaDdlProfiles.Config parityConfig(final RankDualNcrnaWindows.Config saved,
			final String left, final String right, final String suffix){
		assert(saved.models==2 && saved.groups==3) : "The parity witness retains its original two-model/three-input metadata";
		final VerifyNcrnaDdlProfiles.Config c=new VerifyNcrnaDdlProfiles.Config();
		c.profile512=left;
		c.profile1024=right;
		c.panels=saved.panels;
		c.selection=saved.selection;
		c.references=saved.references;
		c.ranks=saved.output+".loci.tsv.gz";
		c.output=saved.output+suffix;
		c.models=2;
		c.groups=3;
		return c;
	}

	private static void parityFault(final RankDualNcrnaWindows.Config saved, final String left,
			final String right, final String kind, final String reason){
		final VerifyNcrnaDdlProfiles.Config c=parityConfig(saved, left, right, "_parity_reject_"+kind+".tsv");
		final String path=saved.output+"_parity_fault_"+kind+".tsv";
		final LineParser1 p=new LineParser1('\t');
		final ByteBuilder b=new ByteBuilder();
		int selected=0;
		try(Input input=new Input(c.ranks); Output output=new Output(path)){
			output.writer.println(input.reader.nextLine());
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				if(!p.termEquals("partition_plain7_masked10_512", 0)){
					output.writer.println(row);
					continue;
				}
				final int index=selected++;
				if(kind.equals("missing") && index==5){continue;}
				if((kind.equals("rank") || kind.equals("input")) && index==0){
					b.clear();
					for(int field=0; field<10; field++){
						if(field>0){b.tab();}
						if(kind.equals("rank") && field==8){b.append(p.parseInt(8)+1);}
						else if(kind.equals("input") && field==3){b.append("wrong_input");}
						else{
							p.setBounds(field);
							b.append(row, p.a(), p.b()-p.a());
						}
					}
					output.writer.print(b.nl());
				}else{output.writer.println(row);}
				if(kind.equals("duplicate") && index==5){output.writer.println(row);}
			}
		}
		check(selected==6, "Each parity fault must start from the complete six-parent selected-method view");
		c.ranks=path;
		rejectParity(c, reason);
	}

	private static void rejectParity(final VerifyNcrnaDdlProfiles.Config c, final String reason){
		try{VerifyNcrnaDdlProfiles.run(c);}
		catch(IllegalArgumentException failure){
			check(failure.getMessage().contains(reason) && !new File(c.output).exists(),
				"Parity must reject before publishing a success table, at the intended mismatch: "+failure.getMessage());
			return;
		}
		throw new AssertionError("Accepted malformed parity input: "+reason);
	}

	/** Two raw masked words deliberately share BOTH the bucket and encoded key.
	 * The independently constructed short witness has exactly one collision bin;
	 * it therefore detects counting compressed keys instead of distinct words. */
	private static void consensusCollisions(){
		final IntHashMap firstWords=new IntHashMap(4096);
		for(int raw=0; raw<(1<<18); raw++){
			final int word=((raw>>>8)<<10)|(raw&255);
			final long hash=Tools.hash64shift(word^12345L);
			final int bucket=(int)(hash&511);
			final int pair=(bucket<<16)|NcrnaDdlQueryKeys.encodeKey(hash);
			final int first=firstWords.get(pair);
			if(first<0){
				firstWords.put(pair, word);
				continue;
			}
			check(first!=word, "The collision witness requires two distinct raw masked words");
			final String a=dna(first), b=dna(word), sequence=a+"N"+b;
			final IntHashSet words=new IntHashSet(16);
			final int[] bins=new int[512];
			check(RankCompactNcrnaWindows.countConsensusWords(sequence.getBytes(StandardCharsets.US_ASCII), words, bins)==2
				&& bins[bucket]==2 && RankCompactNcrnaWindows.countCollisionBins(bins)==1,
				"Two distinct raw words sharing a compressed key must produce exactly one reported collision bin");
			for(int i=0; i<bins.length; i++){
				check(bins[i]==(i==bucket ? 2 : 0), "The two-word witness must occupy exactly its one known collision bin");
			}
			check(RankCompactNcrnaWindows.countConsensusWords((sequence+"N"+a).getBytes(StandardCharsets.US_ASCII), words, bins)==2
				&& bins[bucket]==2 && RankCompactNcrnaWindows.countCollisionBins(bins)==1,
				"Repeating a consensus word cannot add a distinct word or increase reported collision bins");
			final String maskedN=a.substring(0, 5)+"N"+a.substring(6);
			check(RankCompactNcrnaWindows.countConsensusWords(maskedN.getBytes(StandardCharsets.US_ASCII), words, bins)==0
				&& RankCompactNcrnaWindows.countCollisionBins(bins)==0,
				"N at the masked position must invalidate the complete ten-base span");
			for(final int count:bins){check(count==0, "An invalid replacement must clear every prior consensus-word bin");}
			System.err.println("NCRNA_CONSENSUS_COLLISION_WITNESS\t"+sequence+"\tbucket="+bucket+"\traw_words=2\tcollision_bins=1\tcompressed_keys=1");
			return;
		}
		throw new AssertionError("The finite raw-word domain must provide the declared compressed-collision witness");
	}
	private static String dna(final int word){
		final char[] bases=new char[10];
		for(int i=0; i<bases.length; i++){bases[9-i]="ACTG".charAt((word>>>(2*i))&3);}
		return new String(bases);
	}

	private static int rows(final String path){
		int count=0;
		try(Input in=new Input(path)){while(in.reader.nextLine()!=null){count++;}}
		return count;
	}
	private static void equal(final String a, final String b){
		try(Input left=new Input(a); Input right=new Input(b)){
			for(byte[] row; (row=left.reader.nextLine())!=null;){
				check(Arrays.equals(row, right.reader.nextLine()), "The actual compact CLI changed a report row, including UTF-8 fields");
			}
			check(right.reader.nextLine()==null, "CLI cannot append extra compact rows");
		}
	}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
	private static final String[] SUFFIXES={".methods.tsv", ".profile_models.tsv", ".loci.tsv.gz", ".summary.tsv", ".last_losses.tsv", ".inputs.tsv"};
}
