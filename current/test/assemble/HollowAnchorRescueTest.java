package test.assemble;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Random;

import assemble.Tadpole2;
import assemble.HollowAnchorRescue;
import assemble.KmerWindowDepth;
import dna.AminoAcid;
import parse.PreParser;
import shared.Timer;
import stream.Read;
import structures.ByteBuilder;
import ukmer.AbstractKmerTableU;
import ukmer.Kmer;
import ukmer.KmerTableSetU;

/** Tests real extension transactions using controlled in-memory count graphs.
 * No read fixtures are generated; an existing tiny input initializes each table.
 * @author Brian Bushnell, Fischl */
public final class HollowAnchorRescueTest {

	/** Covers word boundaries, both orientations, safe abstention and byte-exact rollback. */
	public static void main(String[] args){
		final PreParser pp=new PreParser(args, HollowAnchorRescueTest.class, false);
		args=pp.args;
		String input=null;
		for(String arg : args){
			final String[] split=arg.split("=", 2);
			if(split.length==2 && split[0].equalsIgnoreCase("in")){input=split[1];}
			else{throw new IllegalArgumentException("Expected in=<existing small fixture>, not "+arg);}
		}
		if(input==null){throw new IllegalArgumentException("An existing input fixture is required.");}
		final Random random=new Random(42124);
		for(int k : new int[]{65, 96, 124, 128, 192}){
			final Tadpole2 baseline=make(input, k, false), rescued=make(input, k, true);
			baseline.loadKmers(new Timer());
			rescued.loadKmers(new Timer());
			for(int stop : new int[]{0, 4, 12}){
				for(int strand=0; strand<2; strand++){
					for(int kind=0; kind<8; kind++){
						exercise(baseline, rescued, random, k, stop, strand, kind);
					}
				}
			}
			baseline.clearData();
			rescued.clearData();
		}
		System.out.println("HOLLOW_ANCHOR_RESCUE_TEST_PASS checks="+checks);
		final int before=checks;
		exerciseSeedGuard(input);
		System.out.println("HOLLOW_SEED_GUARD_INTEGRATION_PASS checks="+(checks-before));
		final int beforePaths=checks;
		exerciseSeedGuardPaths(input);
		System.out.println("HOLLOW_SEED_GUARD_PATH_PASS checks="+(checks-beforePaths));
		final int beforeRelative=checks;
		exerciseRelativeGuard(input);
		System.out.println("HOLLOW_SEED_GUARD_RELATIVE_PASS checks="+(checks-beforeRelative));
		final int beforeProposed=checks;
		exerciseProposedGuard(input);
		System.out.println("HOLLOW_SEED_GUARD_PROPOSED_PASS checks="+(checks-beforeProposed));
		final int beforeDelayed=checks;
		exerciseDelayedGuard(input);
		System.out.println("HOLLOW_DELAYED_ROLLBACK_TEST_PASS checks="+(checks-beforeDelayed));
		final int beforePartial=checks;
		exercisePartialGuard(input);
		System.out.println("HOLLOW_PARTIAL_SEED_GUARD_TEST_PASS checks="+(checks-beforePartial));
		final int beforeCompetitor=checks;
		exerciseCompetitorGuard(input);
		System.out.println("HOLLOW_RESCUE_COMPETITOR_TEST_PASS checks="+(checks-beforeCompetitor));
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fixtures           ----------------*/
	/*--------------------------------------------------------------*/

	/** Explicit CLI settings keep ordinary extension identical between paired engines. */
	private static Tadpole2 make(final String input, final int k, final boolean rescue){
		return make(input, k, rescue, false);
	}

	/** Enables the delayed experiment through the public parser, independently of recovery. */
	private static Tadpole2 make(final String input, final int k, final boolean rescue, final boolean delayed){
		return make(input, k, rescue, delayed, false);
	}

	/** The competitor option uses no supplemental evidence table. */
	private static Tadpole2 make(final String input, final int k, final boolean rescue,
			final boolean delayed, final boolean competitor){
		assert(k>64) : "The left ingress must be in the ignored interior before the next shift.";
		return new Tadpole2(new String[]{"in="+input, "k="+k, "mode=extend", "hollow=t",
				"hollowrescue="+rescue, "hollowbranchrollback="+delayed, "hollowrescuecompetitor="+competitor,
				"t=1", "ways=7", "initialsize=1009", "minprob=0",
				"prefilter=f", "pop=f", "sweeplen=0", "ordered=t", "interleaved=f",
				"extendrollback=0"}, false);
	}

	/** Tests the inclusive support floor, row exclusion, old-mode parity and trial rollback. */
	private static void exerciseCompetitorGuard(final String input){
		final Random random=new Random(421134);
		for(int k : new int[]{96, 128, 192}){
			final Tadpole2 disabled=make(input, k, true), enabled=make(input, k, true, false, true);
			disabled.loadKmers(new Timer()); enabled.loadKmers(new Timer());
			for(int stop : new int[]{0, 4, 12}){
				for(int strand=0; strand<2; strand++){
					for(int depth : new int[]{3, 4}){
						for(boolean sameIngress : new boolean[]{false, true}){
							final int seed=250, checkpoint=seed+stop, position=checkpoint-k+32;
							final byte[] truth=new byte[seed+EXTENSION+1];
							for(int i=0; i<truth.length; i++){truth[i]=AminoAcid.numberToBase[random.nextInt(4)];}
							if(strand!=0){AminoAcid.reverseComplementBasesInPlace(truth);}
							populate(disabled, truth, seed, 100); populate(enabled, truth, seed, 100);
							final byte[] competitor=truth.clone();
							competitor[sameIngress ? checkpoint : position]=alternate(competitor[sameIngress ? checkpoint : position], 2);
							setWindow(disabled, competitor, checkpoint+1, depth);
							setWindow(enabled, competitor, checkpoint+1, depth);
							final byte[] observed=Arrays.copyOf(truth, seed);
							observed[position]=alternate(observed[position], 1);
							final byte[] stopped=Arrays.copyOf(truth, checkpoint), full=Arrays.copyOf(truth, seed+EXTENSION);
							System.arraycopy(observed, 0, stopped, 0, seed);
							final HollowAnchorRescue old=new HollowAnchorRescue(enabled.tables);
							final HollowAnchorRescue explicitOff=new HollowAnchorRescue(enabled.tables, false);
							final HollowAnchorRescue guarded=new HollowAnchorRescue(enabled.tables, true);
							final int replacement=AminoAcid.baseToNumber[truth[position]];
							check(old.propose(stopped, checkpoint, 2)==replacement &&
									explicitOff.propose(stopped, checkpoint, 2)==replacement,
									"Default and explicit-off proposals must retain the old fourfold winner.");
							final boolean veto=!sameIngress && depth==4;
							check(guarded.propose(stopped, checkpoint, 2)==(veto ? -1 : replacement) && guarded.collapsed,
									"Only a different ingress at support>=4 may trigger the additional veto.");
							// A same-ingress competitor can independently stop ordinary extension;
							// its row-exclusion contract is tested directly above, not conflated here.
							if(sameIngress){continue;}
							final ByteBuilder control=prefix(observed);
							check(disabled.extendToRight2(control, null, new int[4], EXTENSION, true, new Kmer(k, true))==EXTENSION &&
									Arrays.equals(control.toBytes(), full), "The off-mode extension must still complete the rescue.");
							final byte[] quality=new byte[seed]; Arrays.fill(quality, (byte)30);
							final Read read=new Read(observed.clone(), quality.clone(), "competitor", checks);
							check(enabled.extendRead(read, new ByteBuilder(), null, new int[4], EXTENSION, new Kmer(k, true))==(veto ? stop : EXTENSION),
									"Public competitor flag must veto before a seed edit or trial extension.");
							check(Arrays.equals(read.bases, veto ? stopped : full), "Competitor veto changed retained sequence bytes.");
							if(!veto){quality[position]=0;}
							check(read.quality.length==read.length() && Arrays.equals(Arrays.copyOf(read.quality, seed), quality),
									"Vetoed corrections must preserve original qualities.");
							if(!veto){
								weakenTail(enabled, truth, 2);
								final ByteBuilder trial=prefix(observed);
								check(enabled.extendToRight2(trial, null, new int[4], EXTENSION, true, new Kmer(k, true))==stop &&
										Arrays.equals(trial.toBytes(), stopped), "Competitor screening must not weaken existing collapse rollback.");
							}
						}
					}
				}
			}
			disabled.clearData(); enabled.clearData();
		}
		boolean rejected=false;
		try{make(input, 128, false, false, true);}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "A competitor veto without hollowrescue must fail at argument validation.");
	}

	/** Sets one complete Hollow window without overwriting shared prefix context. */
	private static void setWindow(final Tadpole2 engine, final byte[] path, final int end, final int count){
		assert(end>=engine.k() && end<=path.length && count>0) : "A controlled competitor requires one full positive-count window.";
		final Kmer word=new Kmer(engine.k(), true);
		for(int i=end-engine.k(); i<end; i++){word.addRight(path[i]);}
		final AbstractKmerTableU[] tables=engine.tables.tables();
		tables[word.mod(tables.length)].set(word, count);
	}

	/** Tests full absolute protection, partial fourfold support and unchanged off mode. */
	private static void exercisePartialGuard(final String input){
		final Random random=new Random(421124);
		for(int k : new int[]{96, 128, 192}){
			final String[] base=guardArgs(input, k), args=Arrays.copyOf(base, base.length+2);
			args[base.length]="hollowseedguardpartial=t";
			args[base.length+1]="hollowseedguardproposed=t";
			final GuardedEngine enabled=new GuardedEngine(args);
			final String[] offArgs=args.clone(); offArgs[base.length]="hollowseedguardpartial=f";
			final GuardedEngine disabled=new GuardedEngine(offArgs);
			enabled.loadKmers(new Timer()); disabled.loadKmers(new Timer());
			for(int stop : new int[]{0, 4, 12}){
				for(int strand=0; strand<2; strand++){
					for(int weak : new int[]{0, 1}){
						final int seed=250, position=seed-k+32+stop;
						final byte[] truth=new byte[seed+EXTENSION+1];
						for(int i=0; i<truth.length; i++){truth[i]=AminoAcid.numberToBase[random.nextInt(4)];}
						if(strand!=0){AminoAcid.reverseComplementBasesInPlace(truth);}
						populate(enabled, truth, seed, 30); populate(disabled, truth, seed, 30);
						final byte[] observed=Arrays.copyOf(truth, seed);
						observed[position]=alternate(observed[position], 1);
						final byte[] full=Arrays.copyOf(truth, seed+EXTENSION), stopped=Arrays.copyOf(truth, seed+stop);
						System.arraycopy(observed, 0, stopped, 0, seed);
						if(weak>0){
							seedSupport(enabled.support(), observed, seed, position, weak);
							seedSupport(disabled.support(), observed, seed, position, weak);
						}
						seedSupport(enabled.support(), truth, seed, position, 1);
						seedSupport(disabled.support(), truth, seed, position, 1);
						ByteBuilder result=prefix(observed);
						check(enabled.extendToRight2(result, null, new int[4], EXTENSION, true, new Kmer(k, true))==EXTENSION &&
								Arrays.equals(result.toBytes(), full), "Absent or depth-one originals must not impose blanket solid support.");
						seedOneWindow(enabled.support(), observed, seed, position, 2);
						seedOneWindow(disabled.support(), observed, seed, position, 2);
						for(int proposed : new int[]{7, 8}){
							seedSupport(enabled.support(), truth, seed, position, proposed);
							seedSupport(disabled.support(), truth, seed, position, proposed);
							final ByteBuilder old=prefix(observed);
							check(disabled.extendToRight2(old, null, new int[4], EXTENSION, true, new Kmer(k, true))==EXTENSION &&
									Arrays.equals(old.toBytes(), full), "Disabled partial protection must retain existing recovery.");
							final byte[] quality=new byte[seed]; Arrays.fill(quality, (byte)30);
							final Read read=new Read(observed.clone(), quality.clone(), "partial", checks);
							final int expected=proposed==8 ? EXTENSION : stop;
							check(enabled.extendRead(read, new ByteBuilder(), null, new int[4], EXTENSION, new Kmer(k, true))==expected,
									"Partial protection must use inclusive fourfold support, including an original minimum of zero.");
							check(Arrays.equals(read.bases, proposed==8 ? full : stopped), "Partial veto must preserve the seed and stopping point.");
							if(proposed==8){quality[position]=0;}
							check(read.quality.length==read.length() && Arrays.equals(Arrays.copyOf(read.quality, seed), quality),
									"Only a committed correction may alter original qualities.");
						}
						weakenTail(enabled, truth, 2);
						result=prefix(observed);
						check(enabled.extendToRight2(result, null, new int[4], EXTENSION, true, new Kmer(k, true))==stop &&
								Arrays.equals(result.toBytes(), stopped), "Partial support must not bypass trial rollback.");
						populate(enabled, truth, seed, 30);
						seedSupport(enabled.support(), observed, seed, position, 2);
						seedSupport(enabled.support(), truth, seed, position, 400);
						result=prefix(observed);
						check(enabled.extendToRight2(result, null, new int[4], EXTENSION, true, new Kmer(k, true))==stop &&
								Arrays.equals(result.toBytes(), stopped), "A repeat-majority proposal must never override full original support.");
					}
				}
			}
			enabled.clearData(); disabled.clearData();
		}
		final String[] base=guardArgs(input, 128), invalid=Arrays.copyOf(base, base.length+2);
		invalid[base.length]="hollowseedguardpartial=t";
		invalid[base.length+1]="hollowseedguardrelative=t"; expectGuardRejected(invalid);
		invalid[base.length+1]="hollowseedguard=null"; expectGuardRejected(invalid);
		invalid[base.length]="hollowseedguardpartial=f";
		final GuardedEngine disabled=new GuardedEngine(invalid);
		check(disabled.support()==null, "Disabled partial protection must not require an evidence table.");
		disabled.clearData();
	}

	/** Supports one original window while leaving other overlapping windows weak. */
	private static void seedOneWindow(final KmerTableSetU tables, final byte[] bases,
			final int length, final int position, final int count){
		final int first=Math.max(0, position-95), last=Math.min(position, length-96);
		assert(tables.kbig==96 && first<last && count>0) :
				"The partial fixture needs multiple seed-contained windows; zero cannot be written as a table deletion.";
		final Kmer word=new Kmer(96, false);
		for(int i=first; i<first+96; i++){word.addRight(bases[i]);}
		tables.tables()[word.mod(tables.ways)].set(word, count);
	}

	/** Exercises completion, pre-rescue rollback and rollback of an already edited trial.
	 * All sequences are in-memory count-graph fixtures, not new synthetic read files. */
	private static void exerciseDelayedGuard(final String input){
		final Random random=new Random(421120);
		for(int k : new int[]{96, 128, 192}){
			final Tadpole2 plain=make(input, k, false, true);
			final Tadpole2 delayed=make(input, k, true, true);
			final Tadpole2 control=make(input, k, true, false);
			plain.loadKmers(new Timer()); delayed.loadKmers(new Timer()); control.loadKmers(new Timer());
			for(int span : new int[]{k-1, k+1}){
				for(int strand=0; strand<2; strand++){
					for(int kind=0; kind<5; kind++){
						// 0 complete/no edit;1 stall/no edit;2 warning before rescue;
						// 3 warning during a stalled trial;4 complete trial with warning.
						final int seed=250, stop=4, position=seed-k+32+stop;
						final byte[] truth=new byte[seed+EXTENSION+1];
						for(int i=0; i<truth.length; i++){truth[i]=AminoAcid.numberToBase[random.nextInt(4)];}
						if(strand!=0){AminoAcid.reverseComplementBasesInPlace(truth);}
						final byte[] path=kind==1 ? Arrays.copyOf(truth, seed+15) : truth;
						final Tadpole2 engine=kind<2 ? plain : delayed;
						populate(engine, path, seed, 30); populate(control, path, seed, 30);
						final byte[] observed=Arrays.copyOf(truth, seed);
						if(kind>=2){observed[position]=alternate(observed[position], 1);}
						final byte[] originalPath=truth.clone();
						System.arraycopy(observed, 0, originalPath, 0, seed);
						final int warning=kind<3 ? 2 : 10;
						warnAt(engine, kind<3 ? originalPath : truth, seed+warning, span);
						if(kind==3){
							branchAt(engine, truth, seed+15); branchAt(control, truth, seed+15);
						}
						final int expected=kind==0 || kind==4 ? EXTENSION : (kind==3 ? stop : warning);
						final byte[] expectedBases=Arrays.copyOf(truth, seed+expected);
						if(kind!=4){System.arraycopy(observed, 0, expectedBases, 0, seed);}
						for(boolean junction : new boolean[]{false, true}){
							final ByteBuilder result=prefix(observed);
							check(engine.extendToRight2(result, null, new int[4], EXTENSION, junction, new Kmer(k, true))==expected,
									"Delayed rollback yield differs: K="+k+", span="+span+", kind="+kind);
							check(Arrays.equals(result.toBytes(), expectedBases), "Delayed rollback changed seed or retained suffix bytes.");
							if(kind==3){
								final ByteBuilder old=prefix(observed);
								final int added=control.extendToRight2(old, null, new int[4], EXTENSION, junction, new Kmer(k, true));
								check(added>stop && added<EXTENSION && old.array[position]==truth[position],
										"The trial control must retain its edit: existing collapse rollback must not explain this result.");
							}
						}
						final byte[] quality=new byte[seed]; Arrays.fill(quality, (byte)30);
						final Read read=new Read(observed.clone(), quality.clone(), "delayed", checks);
						check(engine.extendRead(read, new ByteBuilder(), null, new int[4], EXTENSION, new Kmer(k, true))==expected,
								"Read and builder delayed checkpoints must agree.");
						check(Arrays.equals(read.bases, expectedBases) && read.quality.length==read.length(), "Read reconstruction lost retained bytes or quality length.");
						if(kind==4){quality[position]=0;}
						check(Arrays.equals(Arrays.copyOf(read.quality, seed), quality), "Only a retained seed correction may receive Q0.");
						// Immediate initial-count failure must not reuse a preceding checkpoint.
						final byte[] absent=new byte[seed];
						for(int i=0; i<seed; i++){absent[i]=AminoAcid.numberToBase[random.nextInt(4)];}
						final ByteBuilder empty=prefix(absent);
						check(engine.extendToRight2(empty, null, new int[4], EXTENSION, true, new Kmer(k, true))==0 &&
								Arrays.equals(empty.toBytes(), absent), "Early return leaked a previous read checkpoint.");
					}
				}
			}
			plain.clearData(); delayed.clearData(); control.clearData();
		}
		final String[] invalid=guardArgs(input, 128);
		final String[] both=Arrays.copyOf(invalid, invalid.length+2);
		both[invalid.length]="hollowbranchstop=t"; both[invalid.length+1]="hollowbranchrollback=t";
		expectGuardRejected(both);
	}

	/** Installs a hard alternate-span contradiction without changing nominal path counts. */
	private static void warnAt(final Tadpole2 engine, final byte[] path, final int end, final int span){
		assert(end>=span && end<path.length) : "A warning fixture needs observed context and a proposed next base.";
		final AbstractKmerTableU[] tables=engine.tables.tables();
		final Kmer word=new Kmer(engine.k(), true), next=new Kmer(engine.k(), true);
		KmerWindowDepth.setHollowQuery(word, path, end-span, span);
		tables[word.mod(tables.length)].set(word, 10);
		final int chosen=AminoAcid.baseToNumber[path[end]];
		for(int base=0; base<4; base++){
			next.setFrom(word); next.addRightNumeric(base);
			tables[next.mod(tables.length)].set(next, base==((chosen+1)&3) ? 30 : 1);
		}
	}

	/** Ends a trial with a supported ordinary branch, not a second low-count collapse. */
	private static void branchAt(final Tadpole2 engine, final byte[] path, final int end){
		assert(end>=engine.k() && end<path.length) : "The branch must lie inside the trial suffix.";
		final Kmer word=new Kmer(engine.k(), true);
		for(int i=end-engine.k(); i<end; i++){word.addRight(path[i]);}
		word.addRightNumeric((AminoAcid.baseToNumber[path[end]]+1)&3);
		final AbstractKmerTableU[] tables=engine.tables.tables();
		tables[word.mod(tables.length)].set(word, 30);
	}

	/** Exercises the actual recovery wrapper with and without an independently populated veto table. */
	private static void exerciseSeedGuard(final String input){
		final Random random=new Random(421053);
		for(int k : new int[]{96, 128, 192}){
			final Tadpole2 unguarded=make(input, k, true);
			final GuardedEngine guarded=new GuardedEngine(guardArgs(input, k));
			unguarded.loadKmers(new Timer());
			guarded.loadKmers(new Timer());
			for(int stop : new int[]{0, 4, 12}){
				for(int strand=0; strand<2; strand++){
					final int seed=250, position=seed-k+32+stop;
					final byte[] truth=new byte[seed+EXTENSION+1];
					for(int i=0; i<truth.length; i++){truth[i]=AminoAcid.numberToBase[random.nextInt(4)];}
					if(strand!=0){AminoAcid.reverseComplementBasesInPlace(truth);}
					populate(unguarded, truth, seed, 30);
					populate(guarded, truth, seed, 30);
					final byte[] observed=Arrays.copyOf(truth, seed);
					observed[position]=alternate(observed[position], 1);
					final byte[] full=Arrays.copyOf(truth, seed+EXTENSION);
					ByteBuilder result=prefix(observed);
					check(unguarded.extendToRight2(result, null, new int[4], EXTENSION, true, new Kmer(k, true))==EXTENSION,
							"The unguarded fixture must expose the recoverable proposal.");
					check(Arrays.equals(result.toBytes(), full), "Unguarded recovery differs from the seeded path.");
					result=prefix(observed);
					check(guarded.extendToRight2(result, null, new int[4], EXTENSION, true, new Kmer(k, true))==EXTENSION,
							"Absent solid support must not impose a blanket solid requirement.");
					check(Arrays.equals(result.toBytes(), full), "Unsupported original seed must retain native recovery.");
					final KmerTableSetU support=guarded.support();
					final Kmer word=new Kmer(96, false);
					final int first=Math.max(0, position-95), last=Math.min(position, seed-96);
					for(int start=first; start<=last; start++){
						word.clear();
						for(int i=start; i<start+96; i++){word.addRight(observed[i]);}
						support.tables()[word.mod(support.ways)].set(word, 2);
					}
					final byte[] qualities=new byte[seed];
					Arrays.fill(qualities, (byte)30);
					Read read=new Read(observed.clone(), qualities.clone(), "guarded", checks);
					check(guarded.extendRead(read, new ByteBuilder(), null, new int[4], EXTENSION, new Kmer(k, true))==stop,
							"Supported original seed must veto before trial extension.");
					final byte[] stopped=Arrays.copyOf(truth, seed+stop);
					System.arraycopy(observed, 0, stopped, 0, seed);
					check(Arrays.equals(read.bases, stopped), "Seed veto must preserve the original stopping endpoint and seed bytes.");
					check(Arrays.equals(Arrays.copyOf(read.quality, seed), qualities), "A veto must not lower an unedited seed quality.");
					// word still holds the last overlapping original window. This
					// one weak count must release the veto, even after its prior use.
					support.tables()[word.mod(support.ways)].set(word, 1);
					read=new Read(observed.clone(), qualities.clone(), "weak_last", checks);
					check(guarded.extendRead(read, new ByteBuilder(), null, new int[4], EXTENSION, new Kmer(k, true))==EXTENSION,
							"One weak original window must allow the ordinary recovery transaction.");
					check(Arrays.equals(read.bases, full), "Released recovery changed the expected path.");
					check(read.quality[position]==0, "An accepted correction still receives the existing Q0 annotation.");
					check(support.getCount(word)==1, "Guard queries must not modify the independent evidence table.");
				}
			}
			final KmerTableSetU owned=guarded.support();
			unguarded.clearData(); guarded.clearData();
			check(owned.tables()==null, "clearData must release the independent seed-support partitions too.");
		}
		final java.io.File fixture=new java.io.File(input).getAbsoluteFile();
		final String alias=new java.io.File(fixture.getParentFile(), "./"+fixture.getName()).getPath();
		for(String invalid : new String[]{"hollowrescue=f", "el=1", "extendrollback=3", "hollowseedguard", "hollowseedguard=stdin",
				"outstream="+input, "out="+alias, "--outstream="+input, "-out="+alias}){
			final String[] valid=guardArgs(input, 128), args=Arrays.copyOf(valid, valid.length+1);
			args[valid.length]=invalid;
			expectGuardRejected(args);
		}
		final String[] valid=guardArgs(input, 128), hyphenated=Arrays.copyOf(valid, valid.length+1);
		for(int i=0; i<valid.length; i++){
			if(hyphenated[i].startsWith("hollowseedguard=")){hyphenated[i]="--"+hyphenated[i];}
		}
		hyphenated[valid.length]="--outstream="+input;
		expectGuardRejected(hyphenated);
		final String[] disabled=Arrays.copyOf(valid, valid.length+1);
		disabled[valid.length]="hollowseedguard=null";
		final GuardedEngine off=new GuardedEngine(disabled);
		check(off.support()==null, "Literal null must disable protection after PreParser normalization.");
		off.clearData();
	}

	/** Rejects unsupported settings and evidence/output aliases before creating output. */
	private static void expectGuardRejected(final String[] args){
		boolean rejected=false;
		try{new GuardedEngine(args);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Unsupported guard configuration must fail: "+Arrays.toString(args));
	}

	/** Tests CLI wiring, the fourfold boundary and rollback after releasing a veto. */
	private static void exerciseRelativeGuard(final String input){
		final Random random=new Random(421061);
		for(int k : new int[]{96, 128, 192}){
			final String[] base=guardArgs(input, k), args=Arrays.copyOf(base, base.length+1);
			args[base.length]="hollowseedguardrelative=t";
			final GuardedEngine engine=new GuardedEngine(args);
			engine.loadKmers(new Timer());
			for(int stop : new int[]{0, 4, 12}){
				for(int strand=0; strand<2; strand++){
					final int seed=250, position=seed-k+32+stop;
					final byte[] truth=new byte[seed+EXTENSION+1];
					for(int i=0; i<truth.length; i++){truth[i]=AminoAcid.numberToBase[random.nextInt(4)];}
					if(strand!=0){AminoAcid.reverseComplementBasesInPlace(truth);}
					populate(engine, truth, seed, 30);
					final byte[] observed=Arrays.copyOf(truth, seed);
					observed[position]=alternate(observed[position], 1);
					seedSupport(engine.support(), observed, seed, position, 2);
					for(int proposedDepth : new int[]{7, 8}){
						seedSupport(engine.support(), truth, seed, position, proposedDepth);
						final byte[] quality=new byte[seed];
						Arrays.fill(quality, (byte)30);
						final Read read=new Read(observed.clone(), quality, "relative", checks);
						final int expected=proposedDepth==8 ? EXTENSION : stop;
						check(engine.extendRead(read, new ByteBuilder(), null, new int[4], EXTENSION, new Kmer(k, true))==expected,
								"Live relative veto must use the inclusive fourfold support boundary.");
						final byte[] expectedBases=Arrays.copyOf(truth, seed+expected);
						if(proposedDepth==7){System.arraycopy(observed, 0, expectedBases, 0, seed);}
						check(Arrays.equals(read.bases, expectedBases), "Relative veto changed a retained seed or extension byte.");
						check(read.quality[position]==(proposedDepth==8 ? 0 : 30), "Only a committed correction may receive Q0.");
					}
					weakenTail(engine, truth, 2);
					final ByteBuilder result=prefix(observed);
					check(engine.extendToRight2(result, null, new int[4], EXTENSION, true, new Kmer(k, true))==stop,
							"Relative seed support must not bypass renewed-collapse rollback.");
					final byte[] stopped=Arrays.copyOf(truth, seed+stop);
					System.arraycopy(observed, 0, stopped, 0, seed);
					check(Arrays.equals(result.toBytes(), stopped), "Rollback must restore the original seed after a released veto.");
				}
			}
			engine.clearData();
		}
		final String[] base=guardArgs(input, 128), invalid=Arrays.copyOf(base, base.length+2);
		invalid[base.length]="hollowseedguardrelative=t";
		invalid[base.length+1]="hollowseedguard=null";
		expectGuardRejected(invalid);
	}

	/** Tests the optional proposed-window veto without changing old recovery behavior.
	 * A single missing proposed word must veto before mutation; depth-one support
	 * must suffice, while original protection and renewed-collapse rollback remain. */
	private static void exerciseProposedGuard(final String input){
		final Random random=new Random(421068);
		for(int k : new int[]{96, 128, 192}){
			final String[] base=guardArgs(input, k), args=Arrays.copyOf(base, base.length+2);
			args[base.length]="hollowseedguardrelative=t";
			args[base.length+1]="hollowseedguardproposed=t";
			final GuardedEngine enabled=new GuardedEngine(args);
			final String[] offArgs=args.clone();
			offArgs[base.length+1]="hollowseedguardproposed=f";
			final GuardedEngine disabled=new GuardedEngine(offArgs);
			enabled.loadKmers(new Timer()); disabled.loadKmers(new Timer());
			for(int stop : new int[]{0, 4, 12}){
				for(int strand=0; strand<2; strand++){
					final int seed=250, position=seed-k+32+stop;
					final byte[] truth=new byte[seed+EXTENSION+1];
					for(int i=0; i<truth.length; i++){truth[i]=AminoAcid.numberToBase[random.nextInt(4)];}
					if(strand!=0){AminoAcid.reverseComplementBasesInPlace(truth);}
					populate(enabled, truth, seed, 30); populate(disabled, truth, seed, 30);
					final byte[] observed=Arrays.copyOf(truth, seed);
					observed[position]=alternate(observed[position], 1);
					seedSupport(enabled.support(), observed, seed, position, 1);
					seedSupport(disabled.support(), observed, seed, position, 1);
					final byte[] full=Arrays.copyOf(truth, seed+EXTENSION), stopped=Arrays.copyOf(truth, seed+stop);
					System.arraycopy(observed, 0, stopped, 0, seed);
					for(int state=0; state<3; state++){
						if(state==1){
							seedSupportExceptLast(enabled.support(), truth, seed, position);
							seedSupportExceptLast(disabled.support(), truth, seed, position);
						}else if(state==2){
							seedSupport(enabled.support(), truth, seed, position, 1);
							seedSupport(disabled.support(), truth, seed, position, 1);
						}
						final ByteBuilder oldResult=prefix(observed);
						check(disabled.extendToRight2(oldResult, null, new int[4], EXTENSION, true, new Kmer(k, true))==EXTENSION,
								"Explicitly disabled proposed support must retain ordinary recovery.");
						check(Arrays.equals(oldResult.toBytes(), full), "Disabled proposed support changed the old corrected bytes.");
						final byte[] quality=new byte[seed]; Arrays.fill(quality, (byte)30);
						final Read read=new Read(observed.clone(), quality, "proposed_support", checks);
						final int expected=state==2 ? EXTENSION : stop;
						check(enabled.extendRead(read, new ByteBuilder(), null, new int[4], EXTENSION, new Kmer(k, true))==expected,
								"The enabled veto requires every proposed window, with depth one sufficient.");
						check(Arrays.equals(read.bases, state==2 ? full : stopped), "A rejected proposal must preserve seed and stopping endpoint.");
						final byte[] expectedQuality=new byte[seed]; Arrays.fill(expectedQuality, (byte)30);
						if(state==2){expectedQuality[position]=0;}
						check(Arrays.equals(Arrays.copyOf(read.quality, seed), expectedQuality), "Only an accepted seed correction may change quality.");
					}
					seedSupport(enabled.support(), observed, seed, position, 2);
					ByteBuilder result=prefix(observed);
					check(enabled.extendToRight2(result, null, new int[4], EXTENSION, true, new Kmer(k, true))==stop,
							"Positive proposed support must not bypass the existing original-allele veto.");
					check(Arrays.equals(result.toBytes(), stopped), "Original protection must keep the uncorrected seed.");
					seedSupport(enabled.support(), truth, seed, position, 8);
					result=prefix(observed);
					check(enabled.extendToRight2(result, null, new int[4], EXTENSION, true, new Kmer(k, true))==EXTENSION,
							"The relative exception and positive proposed support must work together.");
					check(Arrays.equals(result.toBytes(), full), "The combined gates changed an accepted correction.");
					weakenTail(enabled, truth, 2);
					result=prefix(observed);
					check(enabled.extendToRight2(result, null, new int[4], EXTENSION, true, new Kmer(k, true))==stop,
							"Positive proposed support must not bypass renewed-collapse rollback.");
					check(Arrays.equals(result.toBytes(), stopped), "Combined-gate rollback must restore the original seed.");
				}
			}
			enabled.clearData(); disabled.clearData();
		}
		final String[] base=guardArgs(input, 128), invalid=Arrays.copyOf(base, base.length+2);
		invalid[base.length]="hollowseedguardproposed=t";
		invalid[base.length+1]="hollowseedguard=null";
		expectGuardRejected(invalid);
		invalid[base.length]="hollowseedguardproposed=f";
		final GuardedEngine off=new GuardedEngine(invalid);
		check(off.support()==null, "Disabled proposed support must not require or allocate an evidence table.");
		off.clearData();
	}

	/** Leaves the final proposed word absent without breaking a table probe chain. */
	private static void seedSupportExceptLast(final KmerTableSetU tables, final byte[] bases, final int length, final int position){
		assert(tables.kbig==96 && length>=96 && length<=bases.length && position>=0 && position<length) :
				"The missing-word fixture uses the same seed-contained Solid96 coordinates as the live veto.";
		final Kmer word=new Kmer(96, false);
		final int first=Math.max(0, position-95), last=Math.min(position, length-96);
		for(int start=first; start<last; start++){
			word.clear();
			for(int i=start; i<start+96; i++){word.addRight(bases[i]);}
			tables.tables()[word.mod(tables.ways)].set(word, 1);
		}
	}

	/** Installs positive Solid96 counts for every window covering the seed edit. */
	private static void seedSupport(final KmerTableSetU tables, final byte[] bases, final int length,
			final int position, final int count){
		assert(tables.kbig==96 && length>=96 && length<=bases.length && position>=0 && position<length && count>0) :
				"Live fixtures must populate the same original-seed-only Solid96 evidence queried by the guard.";
		final Kmer word=new Kmer(96, false);
		for(int start=Math.max(0, position-95), last=Math.min(position, length-96); start<=last; start++){
			word.clear();
			for(int i=start; i<start+96; i++){word.addRight(bases[i]);}
			tables.tables()[word.mod(tables.ways)].set(word, count);
		}
	}

	/** Tests physical aliases exposed only by mate-name expansion, without risking
	 * the original fixture. Exact temporary paths are removed only after byte parity;
	 * a failed test leaves its small copied fixture available for diagnosis. */
	private static void exerciseSeedGuardPaths(final String input){
		assert(input!=null) : "Path tests copy an existing fixture, never generate synthetic reads.";
		try{
			final Path directory=Files.createTempDirectory("fischl-seed-guard-paths-");
			final Path evidence=Files.copy(Paths.get(input), directory.resolve("evidence.fa"));
			final byte[] before=Files.readAllBytes(evidence);
			for(String output : new String[]{"out", "outd"}){
				for(int mate=1; mate<=2; mate++){
					for(boolean symbolic : new boolean[]{false, true}){
						final Path alias=directory.resolve("alias"+mate+".fa");
						if(symbolic){Files.createSymbolicLink(alias, evidence);}
						else{Files.createLink(alias, evidence);}
						final String[] valid=guardArgs(evidence.toString(), 128);
						final String[] args=Arrays.copyOf(valid, valid.length+1);
						args[valid.length]=output+"="+directory.resolve("alias#.fa");
						boolean rejected=false;
						try{
							final GuardedEngine accepted=new GuardedEngine(args);
							accepted.clearData();
						}catch(IllegalArgumentException expected){
							rejected=expected.getMessage().contains("output would overwrite");
						}
						check(rejected, "Expanded "+output+" mate"+mate+" must reject a physical evidence alias; symbolic="+symbolic);
						check(Arrays.equals(before, Files.readAllBytes(evidence)), "Rejected output must leave evidence bytes unchanged.");
						Files.delete(alias);
					}
				}
			}
			final String[] valid=guardArgs(evidence.toString(), 128);
			final String[] safe=Arrays.copyOf(valid, valid.length+2);
			safe[valid.length]="out="+directory.resolve("safe#.fa");
			safe[valid.length+1]="outd="+directory.resolve("discard#.fa");
			final GuardedEngine accepted=new GuardedEngine(safe);
			accepted.clearData();
			check(Arrays.equals(before, Files.readAllBytes(evidence)), "Distinct mate/discard paths must be accepted without touching evidence.");
			try(java.util.stream.Stream<Path> contents=Files.list(directory)){
				check(contents.count()==1, "Construction must not create read outputs before extension starts.");
			}
			Files.delete(evidence);
			Files.delete(directory);
		}catch(IOException failure){throw new RuntimeException("Could not exercise seed-guard path aliases.", failure);}
	}

	/** Keeps the normal loader and guard loader on existing input, never generated read files. */
	private static String[] guardArgs(final String input, final int k){
		assert(input!=null && k>=96) : "The integration fixture needs an existing input and a full Solid96 seed.";
		return new String[]{"in="+input, "k="+k, "mode=extend", "hollow=t", "hollowrescue=t",
				"hollowseedguard="+input, "el=0", "er=20", "extendrollback=0", "t=1", "ways=7",
				"initialsize=1009", "minprob=0", "prefilter=f", "pop=f", "sweeplen=0", "interleaved=f"};
	}

	/** Test-only access to fill immutable-during-extension counts before each transaction. */
	private static final class GuardedEngine extends Tadpole2 {
		GuardedEngine(final String[] args){super(args, false);}
		KmerTableSetU support(){return hollowSeedGuardTables;}
	}

	/** kind: clean, one recoverable error, repeated collapse, ambiguity, weak context,
	 * terminal depth2, terminal depth3, clean with terminal depth2. */
	private static void exercise(final Tadpole2 baseline, final Tadpole2 rescued,
			final Random random, final int k, final int stop, final int strand, final int kind){
		assert(stop>=0 && stop+5<EXTENSION) : "Both first and recurrent collapses must fit inside the requested extension.";
		// K65 has only one ignored base: later interior positions already belong
		// to its right arm. Use a larger seed, but only test immediate ingress
		// error there; the other geometries test delayed and repeated errors.
		if(k==65 && (stop>0 || kind==2)){return;}
		final int seed=k+1, position=seed-k+32+stop;
		final byte[] truth=new byte[seed+EXTENSION+1];
		for(int i=0; i<truth.length; i++){truth[i]=AminoAcid.numberToBase[random.nextInt(4)];}
		if(strand!=0){AminoAcid.reverseComplementBasesInPlace(truth);}
		final int depth=kind==4 ? 3 : 30;
		populate(baseline, truth, seed, depth);
		populate(rescued, truth, seed, depth);
		if(kind>=5){
			final int tailDepth=kind==6 ? 3 : 2;
			weakenTail(baseline, truth, tailDepth);
			weakenTail(rescued, truth, tailDepth);
		}
		final byte[] observed=Arrays.copyOf(truth, seed);
		if(kind>0 && kind!=7){observed[position]=alternate(observed[position], 1);}
		if(kind==2){observed[position+5]=alternate(observed[position+5], 1);}
		if(kind==3){
			final byte[] decoy=truth.clone();
			decoy[position]=alternate(truth[position], 2);
			populate(baseline, decoy, seed, depth);
			populate(rescued, decoy, seed, depth);
		}
		final int plainExpected=kind==0 || kind==7 ? EXTENSION : stop;
		final int rescueExpected=kind==0 || kind==1 || kind==7 ? EXTENSION : stop;
		for(boolean junction : new boolean[]{false, true}){
			final ByteBuilder plain=prefix(observed), result=prefix(observed);
			check(baseline.extendToRight2(plain, null, new int[4], EXTENSION, junction,
					new Kmer(k, true))==plainExpected, "Baseline stop differs from the controlled count graph.");
			check(rescued.extendToRight2(result, null, new int[4], EXTENSION, junction,
					new Kmer(k, true))==rescueExpected, "Recovery yield/rollback differs from the fixture contract.");
			final byte[] expected=Arrays.copyOf(truth, seed+rescueExpected);
			if(kind!=1){System.arraycopy(observed, 0, expected, 0, seed);}
			check(Arrays.equals(result.toBytes(), expected), "Rollback or abstention changed a retained byte.");
		}
		final byte[] qualities=new byte[seed];
		Arrays.fill(qualities, (byte)30);
		final Read read=new Read(observed.clone(), qualities.clone(), "rescue_fixture", checks);
		check(rescued.extendRead(read, new ByteBuilder(), null, new int[4], EXTENSION,
				new Kmer(k, true))==rescueExpected, "Read and ByteBuilder entry points disagree.");
		check(read.quality.length==read.length(), "Quality length must follow the committed extension only.");
		for(int i=0; i<seed; i++){
			check(read.quality[i]==(kind==1 && i==position ? 0 : qualities[i]),
					"Only the committed substitution may change an original quality.");
		}
		final Read ambiguous=new Read(observed.clone(), null, "undefined", checks);
		ambiguous.bases[seed-2]='N';
		final byte[] unchanged=ambiguous.bases.clone();
		check(rescued.extendRead(ambiguous, new ByteBuilder(), null, new int[4], EXTENSION,
				new Kmer(k, true))==0 && Arrays.equals(ambiguous.bases, unchanged),
				"An undefined terminal window must not enter recovery or reuse previous transaction state.");
	}

	/** Installs full independently constructed windows, including one lookahead step. */
	private static void populate(final Tadpole2 engine, final byte[] path, final int seed, final int count){
		assert(count>0 && seed>=engine.k()) : "Fixture windows must be full and positive; zero is a table sentinel.";
		final AbstractKmerTableU[] tables=engine.tables.tables();
		final Kmer word=new Kmer(engine.k(), true);
		for(int end=seed; end<=path.length; end++){
			word.clear();
			for(int i=end-engine.k(); i<end; i++){word.addRight(path[i]);}
			tables[word.mod(tables.length)].set(word, count);
		}
	}

	/** Makes the final three emitted windows and lookahead weak but traversable.
	 * Ordinary extension must retain these; only an unconfirmed tentative seed
	 * correction should roll back, including at the requested extension limit. */
	private static void weakenTail(final Tadpole2 engine, final byte[] path, final int count){
		assert(count==2 || count==3) : "The regression needs traversable support below the correction confirmation floor4.";
		final AbstractKmerTableU[] tables=engine.tables.tables();
		final Kmer word=new Kmer(engine.k(), true);
		for(int end=path.length-3; end<=path.length; end++){
			word.clear();
			for(int i=end-engine.k(); i<end; i++){word.addRight(path[i]);}
			tables[word.mod(tables.length)].set(word, count);
		}
	}

	/** Changes a base deterministically without introducing undefined sequence. */
	private static byte alternate(final byte base, final int delta){
		assert(delta>0 && delta<4) : "Fixture mutation must change the observed base.";
		return AminoAcid.numberToBase[(AminoAcid.baseToNumber[base]+delta)&3];
	}

	/** Makes an independent mutable seed for each transaction. */
	private static ByteBuilder prefix(final byte[] seed){
		final ByteBuilder builder=new ByteBuilder(seed.length);
		builder.append(seed);
		return builder;
	}

	/** Tests must remain effective even if a caller accidentally disables -ea. */
	private static void check(final boolean pass, final String message){
		if(!pass){throw new AssertionError(message);}
		checks++;
	}

	private static int checks;
	private static final int EXTENSION=20;
}
