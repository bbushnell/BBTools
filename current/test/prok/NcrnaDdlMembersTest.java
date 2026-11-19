package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.Locale;

import cardinality.DynamicDemiLog;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlMembers;
import prok.NcrnaDdlModelProfiles;
import prok.NcrnaDdlProfile;
import prok.NcrnaDdlQueryKeys;
import prok.NcrnaMemberLengthBounds;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Independent literal construction-member joins, including retained cross-labels
 * and checked nonvoting small/orphan records. No learned-resource result.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaDdlMembersTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaDdlMembersTest.class, false);
		final int previous=DynamicDemiLog.exponentBits();
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown fixture parameter: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Supply only out=fresh_fixture_directory and normal runtime options");
			boolean assertions=false;
			assert(assertions=true);
			check(assertions && !Shared.TRIM_READ_DESCRIPTION, "Member fixtures require assertions and complete native read headers");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Fixture output directory must be fresh");
			DynamicDemiLog.setExponent(4);
			run(root.getPath());
			pp.outstream.println("NCRNA_DDL_MEMBERS_TEST_PASS");
		}finally{
			DynamicDemiLog.setExponent(previous);
			Shared.closeStream(pp.outstream);
		}
	}

	private static void run(final String root){
		final NcrnaDdlMembers good=new Fixture(root, "good").load();
		check(good.inputMembers()==5 && good.keptMembers()==3 && good.smallMembers()==1
			&& good.orphanMembers()==1 && good.models()==2, "All five original members retain their independent voting dispositions");
		check(good.members(0).length==2 && good.members(1).length==1,
			"Published construction counts bind profiles to current model order");
		check(Arrays.equals(good.members(0)[0], C.getBytes(StandardCharsets.US_ASCII))
			&& Arrays.equals(good.members(0)[1], A.getBytes(StandardCharsets.US_ASCII))
			&& Arrays.equals(good.members(1)[0], B.getBytes(StandardCharsets.US_ASCII)),
			"Cross-labelled sequence is retained, model groups ignore ledger/assignment order, and members retain FASTA order");
		checkModelProfiles(good);
		forwardProfiles(good);
		consensusProfiles(good);
		final String[] faults={"missing", "duplicate", "description", "undefined", "unknown_model",
			"duplicate_assignment", "missing_assignment", "duplicate_ledger", "wrong_count", "wrong_extrema", "additions", "unknown_status"};
		final String[] reasons={"every routed/assigned", "FASTA full header", "FASTA full header", "undefined base", "absent from the current library",
			"one original builder assignment", "published member summary", "unique complete headers", "census must conserve", "lengths differ from census",
			"without singleton additions", "Unknown builder disposition"};
		check(faults.length==reasons.length, "Each finite defect must name its intended rejection");
		for(int i=0; i<faults.length; i++){
			final Fixture bad=new Fixture(root, faults[i]);
			try{bad.load();}
			catch(IllegalArgumentException failure){
				check(failure.getMessage().contains(reasons[i]), "Wrong rejection for "+faults[i]+": "+failure.getMessage());
				continue;
			}
			throw new AssertionError("Invalid member join was accepted: "+faults[i]);
		}
	}

	/** Literal A/C/G members have one canonical key each. C and G share a
	 * reverse-complement key, so the independent raw-key control determines the
	 * expected one-vote/two-vote buckets without rerunning the profile builder. */
	private static void checkModelProfiles(final NcrnaDdlMembers members){
		for(final int k:new int[]{5, 7, 9, 11, 13, 15, 19}){
			for(int width=32; width<=2048; width<<=1){
				final NcrnaDdlProfile[] profiles=NcrnaDdlModelProfiles.build(members, k, width);
				final int a=homopolymerKey(A, k, width), g=homopolymerKey(B, k, width);
				check(profiles.length==2 && profiles[0].members==2 && profiles[1].members==1,
					"Only the three original kept members contribute, in current model order");
				for(int bucket=0; bucket<width; bucket++){
					final int aKey=(a>>>16)==bucket ? a&65535 : 0;
					final int gKey=(g>>>16)==bucket ? g&65535 : 0;
					final int expected=aKey==0 ? gKey : gKey==0 ? aKey : Math.min(aKey, gKey);
					check(profiles[0].keys()[bucket]==expected && profiles[1].keys()[bucket]==gKey,
						"Literal member votes must match canonical keys: k="+k+", width="+width+", bucket="+bucket);
				}
			}
		}
		final NcrnaDdlProfile[] empty=NcrnaDdlModelProfiles.build(members, 31, 32);
		check(empty[0].members==2 && empty[1].members==1,
			"Members shorter than k retain their membership count while supplying no nonempty keys");
		for(final NcrnaDdlProfile profile:empty){
			for(int bucket=0; bucket<profile.buckets(); bucket++){
				check(profile.keys()[bucket]==0, "No complete k-mer must leave every profile bin empty");
			}
		}
		check(Arrays.equals(members.members(0)[0], C.getBytes(StandardCharsets.US_ASCII)),
			"Native profile construction must borrow retained construction sequences read-only");
	}

	private static int homopolymerKey(final String sequence, final int k, final int width){
		final NcrnaDdlQueryKeys raw=new NcrnaDdlQueryKeys(width, k);
		raw.build(sequence.getBytes(StandardCharsets.US_ASCII));
		final int[] keys=raw.packedKeys();
		check(keys.length==1, "A complete homopolymer has one distinct canonical native bucket/key");
		return keys[0];
	}

	/** Literal packed A/C/G words independently check actual member-profile ingestion. */
	private static void forwardProfiles(final NcrnaDdlMembers members){
		// Native DNA codes are A=0,C=1,G=3. Even k8 masks the right central base.
		final long[][] words={{0x1555L, 0x3fffL, 0x1515L, 0x3f3fL},
			{0x5555L, 0xffffL, 0x5515L, 0xff3fL}};
		for(int k=7; k<=8; k++){
			for(final boolean middle:new boolean[]{false, true}){
				for(final int width:new int[]{256, 2048}){
					final NcrnaDdlProfile[] profiles=NcrnaDdlModelProfiles.build(members, k, width, false, middle);
					final DynamicDemiLog a=DynamicDemiLog.create(width, k, 12345, 0, false);
					final DynamicDemiLog c=DynamicDemiLog.create(width, k, 12345, 0, false);
					final DynamicDemiLog g=DynamicDemiLog.create(width, k, 12345, 0, false);
					a.add(0);
					c.add(words[k-7][middle ? 2 : 0]);
					g.add(words[k-7][middle ? 3 : 1]);
					final NcrnaDdlProfile.Votes expected=new NcrnaDdlProfile.Votes(width);
					expected.add(a.maxArray());
					expected.add(c.maxArray());
					final NcrnaDdlProfile first=expected.finish();
					check(Arrays.equals(profiles[0].keys(), first.keys()) && Arrays.equals(profiles[1].keys(), g.maxArray())
						&& profiles[0].members==2 && profiles[1].members==1,
						"Actual forward member profiles must match literal plain/masked packed A/C/G words");
				}
			}
		}
	}

	/** Literal masked10 homopolymer words check member/consensus ingestion at
	 * all compact widths. Consensus T is absent from the A/C/G member pool.
	 */
	private static void consensusProfiles(final NcrnaDdlMembers members){
		final byte[] consensus="TTTTTTTTTTTTTTTTTTTTTT".getBytes(StandardCharsets.US_ASCII);
		final int[] weights={0, 1, 4, 16, NcrnaDdlProfile.CONSENSUS_WINS};
		for(final int width:new int[]{512, 1024, 2048}){
			final NcrnaDdlProfile[][] profiles=NcrnaDdlModelProfiles.buildConsensus(members,
				new byte[][]{consensus, consensus}, width, weights);
			final char[][] literal=new char[4][];
			final long[] words={0, 0x55455L, 0xffcffL, 0xaa8aaL};
			for(int i=0; i<words.length; i++){
				final DynamicDemiLog sketch=DynamicDemiLog.create(width, 10, 12345, 0, false);
				sketch.add(words[i]);
				literal[i]=sketch.maxArray();
			}
			for(int variant=0; variant<weights.length; variant++){
				for(int model=0; model<2; model++){
					final NcrnaDdlProfile profile=profiles[variant][model];
					final int count=model==0 ? 2 : 1;
					final int weight=weights[variant]<0 ? count+1 : weights[variant];
					check(profile.members==count, "Consensus additions must not change actual model membership");
					for(int bucket=0; bucket<width; bucket++){
						int expected=0, best=0;
						for(final char[] candidate:literal){
							final int key=candidate[bucket];
							if(key==0){continue;}
							final int real=model==0 ? (literal[0][bucket]==key ? 1 : 0)+(literal[1][bucket]==key ? 1 : 0)
								: literal[2][bucket]==key ? 1 : 0;
							final int score=real+(literal[3][bucket]==key ? weight : 0);
							if(score>best || score==best && score>0 && key>expected){
								expected=key;
								best=score;
							}
						}
						check(profile.keys()[bucket]==expected,
							"Masked10 construction must equal literal A/C/G member votes and the declared T consensus bonus");
					}
				}
			}
		}
		check(Arrays.equals(consensus, "TTTTTTTTTTTTTTTTTTTTTT".getBytes(StandardCharsets.US_ASCII)),
			"Consensus construction borrows model sequences read-only");
	}

	private static final class Fixture{
		Fixture(final String root, final String fault){
			dir=root+'/'+fault;
			check(new File(dir).mkdir(), "Each member-join case needs a fresh directory");
			final String header="locus_id\tsource_family\tmember_family\tbuild_group\trouting_rule\trouting_bits\tlength\tinput_header\n";
			final String a="RF00174_a\tRF00174\tRF00174\tRF00174\tACCEPTED_RF00174_PRODUCT\tNA\t20\t"+RAW_A+'\n';
			write("ledger.tsv", header+a
				+"RF01482_c\tRF01482\tB12clade\tRF00174\tHIGHEST_BITS_OPERATIONAL_ONLY\t30\t24\t"+RAW_C+'\n'
				+"RF00174_b\tRF00174\tRF00174\tRF00174\tACCEPTED_RF00174_PRODUCT\tNA\t22\t"+RAW_B+'\n'
				+"RF00174_s\tRF00174\tRF00174\tRF00174\tACCEPTED_RF00174_PRODUCT\tNA\t18\t"+RAW_S+'\n'
				+"RF00174_o\tRF00174\tRF00174\tRF00174\tACCEPTED_RF00174_PRODUCT\tNA\t19\t"+RAW_O+'\n'
				+"RF01689_x\tRF01689\tRF01689\tRF01689\tOWN_SUBTYPE\t40\t20\tRF01689_x source=tid_9_species member_family=RF01689\n"
				+(fault.equals("duplicate_ledger") ? a : ""));
			final String assignedA=assignment(H_A, fault.equals("unknown_model") ? "cobalamin_RF00174_absent" : NAMES[0],
				fault.equals("unknown_status") ? "other" : "kept");
			write("assignments.tsv", "#fullRecordId\tsourceTid\tbuilderGroup\tclusterKey\tmodelLabel\tstatus\n"
				+assignedA+assignment(H_B, NAMES[1], "kept")
				+assignment(H_S, "", "small")+assignment(H_C, NAMES[0], "kept")
				+(fault.equals("missing_assignment") ? "" : assignment(H_O, "", "orphan"))
				+(fault.equals("duplicate_assignment") ? assignedA : ""));
			write("members.fa", fasta(H_B, B)+fasta(H_C, fault.equals("undefined") ? 'N'+C.substring(1) : C)
				+fasta(fault.equals("description") ? H_A+" changed" : H_A, A)+fasta(H_S, S)
				+(fault.equals("missing") ? "" : fasta(H_O, O))+(fault.equals("duplicate") ? fasta(H_A, A) : ""));
			write("models.tsv", NcrnaMemberLengthBounds.HEADER+'\n'
				+"cobalamin_RF00174_fixture0\tRF00174\t22\t"+(fault.equals("wrong_count") ? 3 : 2)+"\t"+(fault.equals("wrong_extrema") ? 19 : 20)+"\t20\t20\t20\t24\t0.5\tMIXED\t1\t0.5\t0\t0\t0\t0\t1\t0.5\n"
				+"cobalamin_RF00174_fixture1\tRF00174\t22\t1\t22\t22\t22\t22\t22\t1\tPURE\t1\t1\t0\t0\t0\t0\t0\t0\n");
			write("summary.tsv", "emitted_family\tinput_members\toriginal_kept\toriginal_small\toriginal_orphan\tadded_memberships\tconstruction_memberships\tmodels\tpure\tmixed\tempty\tpurity_cut\tgroup_min\tgroup_p10\tgroup_median\tgroup_p90\tgroup_max\n"
				+"RF00174\t5\t3\t1\t1\t"+(fault.equals("additions") ? 1 : 0)+"\t3\t2\t1\t1\t0\t0.8\t18\t18\t20\t22\t24\n");
		}
		NcrnaDdlMembers load(){
			return NcrnaDdlMembers.load(dir+"/members.fa", dir+"/ledger.tsv", dir+"/assignments.tsv",
				dir+"/models.tsv", dir+"/summary.tsv", NAMES,
				new byte[][]{B.getBytes(StandardCharsets.US_ASCII), B.getBytes(StandardCharsets.US_ASCII)});
		}
		void write(final String file, final String contents){
			try(Output out=new Output(dir+'/'+file)){out.writer.print(contents.getBytes(StandardCharsets.UTF_8));}
		}
		final String dir;
	}

	private static String assignment(final String header, final String model, final String state){
		return header+"\t1\tunknown\tRF00174_fixture\t"+model+'\t'+state+'\n';
	}
	private static String fasta(final String header, final String bases){return '>'+header+'\n'+bases+'\n';}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private static final String RAW_A="RF00174_a source=tid_1_species status=ANNOTATION_FULL_COPY";
	private static final String RAW_B="RF00174_b source=tid_2_species status=ANNOTATION_FULL_COPY";
	private static final String RAW_C="RF01482_c source=tid_3_species source_family=RF01482 member_family=B12clade";
	private static final String RAW_S="RF00174_s source=tid_4_species status=ANNOTATION_FULL_COPY";
	private static final String RAW_O="RF00174_o source=tid_5_species status=ANNOTATION_FULL_COPY";
	private static final String H_A=RAW_A+" member_family=RF00174 build_group=RF00174";
	private static final String H_B=RAW_B+" member_family=RF00174 build_group=RF00174";
	private static final String H_C=RAW_C+" build_group=RF00174";
	private static final String H_S=RAW_S+" member_family=RF00174 build_group=RF00174";
	private static final String H_O=RAW_O+" member_family=RF00174 build_group=RF00174";
	private static final String A="AAAAAAAAAAAAAAAAAAAA", B="GGGGGGGGGGGGGGGGGGGGGG", C="CCCCCCCCCCCCCCCCCCCCCCCC";
	private static final String S="TTTTTTTTTTTTTTTTTT", O="TTTTTTTTTTTTTTTTTTT";
	private static final String[] NAMES={"cobalamin_RF00174_fixture0 n=2", "cobalamin_RF00174_fixture1 n=1"};
}
