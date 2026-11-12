package prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.FileFormat;
import map.ObjectMap;
import parse.Parse;
import prok.GradeNcrnaReference.Binding;
import prok.NcrnaPolicySearch.Measurement;
import prok.NcrnaTextIO.Output;
import prok.SummarizeNcrnaSweep.Table;
import shared.Tools;
import structures.ByteBuilder;

/** Validates opt-in search coverage from native saved counts, without calling or
 * rereading genomes/GFFs. One driver owns this instance across rounds. The first
 * validated shipped vector remains immutable; equal aggregate TP is insufficient.
 * Input/policy/member rows must follow the native writer's declared order.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaSearchCoverage{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaSearchCoverage(final String[] caller, final String reference, final String[] members_){
		this(caller, reference, members_, false);
	}

	/** Uses the driver's selected logical-work unit for native count conservation. */
	public NcrnaSearchCoverage(final String[] caller, final String reference, final String[] members_, final boolean modelBaseCost){
		metrics=METRICS.clone();
		if(modelBaseCost){metrics[COST]="aligned_model_bases";}
		if(caller==null || reference==null || members_==null){
			throw new IllegalArgumentException("Coverage validation requires caller inputs, references and an explicit member roster");
		}
		members=members_.clone();
		require(members.length<=63, "Coverage roster must fit the qualified-support bit representation");
		bindings=bindInputs(caller, reference, false);
		inputTokens=new byte[bindings.length][];
		for(int i=0; i<bindings.length; i++){inputTokens[i]=token(bindings[i].id);}
		memberTokens=new byte[members.length][];
		for(int m=0; m<members.length; m++){
			require(members[m]!=null && !members[m].isEmpty(), "Member coverage needs named reference members");
			for(int n=0; n<m; n++){require(!members[m].equals(members[n]), "Duplicate coverage member");}
			memberTokens[m]=token(members[m]);
		}
		baseline=new long[Math.multiplyExact(bindings.length, 3)];
		baselineMembers=new long[Math.multiplyExact(bindings.length, members.length*2)];
	}

	/** Exact native caller population/order; off-panel reference payloads stay unopened. */
	static Binding[] bindInputs(final String[] caller, final String reference, final boolean metadataOnly){
		require(caller!=null && reference!=null, "Saved metadata must bind actual caller inputs and references");
		final ArrayList<String> paths=new ArrayList<String>();
		for(final String arg:caller){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("in") || a.equals("infna") || a.equals("fnain") || a.equals("fna") || a.equals("ref")){
				if(metadataOnly){
					if(b==null || b.isEmpty() || b.equalsIgnoreCase("null")){paths.clear();}
					else{
						// One-time config expansion. Exact saved declarations replace live
						// existence/suffix guessing; every path must join the saved reference.
						for(final String path:b.split(",", -1)){paths.add(path);}
					}
				}else{require(Tools.addFiles(b, paths), "Coverage inputs must resolve through the native caller: "+arg);}
			}else if(b==null && (metadataOnly || new File(arg).exists()) && FileFormat.isFastaFile(arg)){paths.add(arg);}
		}
		if(!metadataOnly){Tools.fixExtension(paths);}
		require(!paths.isEmpty(), "Coverage protection needs the actual caller input population");
		final ObjectMap<String, Binding> byPath=new ObjectMap<String, Binding>(String.class, Binding.class);
		for(final Binding binding:GradeNcrnaReference.loadManifest(reference)){byPath.put(binding.fasta, binding);}
		final Binding[] bindings=new Binding[paths.size()];
		for(int i=0; i<bindings.length; i++){
			final Binding binding=byPath.remove(paths.get(i));
			require(binding!=null, "Caller input must bind exactly once to its reference: "+paths.get(i));
			bindings[i]=binding;
		}
		// References outside the caller population deliberately remain unopened.
		return bindings;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Validation         ----------------*/
	/*--------------------------------------------------------------*/

	/** The saved shipped vector detects same-total coverage drift on resume.
	 * It is written once after accounting passes, before the round completion marker. */
	public void bindBaseline(final String path, final boolean resume){
		require(path!=null && baselinePath==null, "Bind one immutable shipped-coverage vector per driver");
		baselinePath=path;
		if(!resume){
			require(!new File(path).exists(), "New search cannot overwrite a shipped-coverage baseline");
			return;
		}
		try(Table table=new Table(path)){
			final int input=table.column("input_id");
			final int[] union=columns(table, BASELINE_COLUMNS);
			final int[] member=columns(table, memberColumns());
			for(int i=0; i<bindings.length; i++){
				require(table.next() && table.matches(inputTokens[i], input), "Saved shipped vector must bind every exact input in order");
				for(int c=0; c<union.length; c++){baseline[i*3+c]=count(table, union[c]);}
				for(int c=0; c<member.length; c++){baselineMembers[i*member.length+c]=count(table, member[c]);}
			}
			require(!table.next(), "Unexpected extra shipped-coverage baseline row");
		}
		bound=true;
	}

	/** Completed-marker publication follows this check, including on resume. */
	public Counts[] read(final String prefix, final String[] policies, final Measurement[] measured){
		require(policies!=null && measured!=null && policies.length==measured.length && policies.length>0,
			"Coverage and objective measurements must share one complete policy roster");
		final byte[][] policyTokens=new byte[policies.length][];
		final Counts[] counts=new Counts[policies.length];
		for(int p=0; p<policies.length; p++){
			policyTokens[p]=token(policies[p]);
			counts[p]=new Counts(0, 0, new long[members.length], new long[members.length]);
		}
		final long[][] totals=new long[policies.length][metrics.length];
		try(Table table=new Table(prefix+".inputs.tsv");
			Table memberTable=members.length==0 ? null : new Table(prefix+".coverage.tsv")){
			final int input=table.column("input_id"), policy=table.column("policy");
			final int[] columns=columns(table, metrics);
			final long[] root=new long[metrics.length], row=new long[metrics.length];
			final long[] rootVector=new long[3], unionTp=new long[policies.length];
			final long[] memberRoot=new long[members.length*2], memberRow=new long[members.length*2];
			for(int i=0; i<bindings.length; i++){
				for(int p=0; p<policies.length; p++){
					require(table.next() && table.matches(inputTokens[i], input) && table.matches(policyTokens[p], policy),
						"Missing, repeated or unexpected input/policy row; expected "+bindings[i].id+" / "+policies[p]);
					for(int c=0; c<columns.length; c++){
						row[c]=count(table, columns[c]);
						totals[p][c]=Math.addExact(totals[p][c], row[c]);
					}
					require(Math.addExact(row[TP], row[FN])==row[SCORED], "Input TP+FN must conserve scored references");
					if(p==0){
						System.arraycopy(row, 0, root, 0, row.length);
						require(row[LOST]==0 && row[GAINED]==0, "Shipped input rows cannot have root-relative changes");
						rootVector[0]=row[TP];
						rootVector[1]=row[SCORED];
						rootVector[2]=row[MASKED];
						checkBaseline(baseline, i*3, rootVector);
					}
					require(row[SCORED]==root[SCORED] && row[MASKED]==root[MASKED] && row[LOST]<=root[TP]
						&& row[GAINED]<=root[SCORED]-root[TP] && Math.addExact(root[TP]-row[LOST], row[GAINED])==row[TP],
						"Root minus lost plus gained must explain each input's TP under unchanged denominators");
					if(root[TP]>0 && row[TP]==0){counts[p].lastLoss++;}
					if(root[TP]==0 && row[TP]>0){counts[p].firstGain++;}
					unionTp[p]=row[TP];
				}
				if(memberTable!=null){readMembers(memberTable, i, policyTokens, counts, unionTp, root[SCORED], memberRoot, memberRow);}
			}
			require(!table.next(), "Unexpected extra input/policy row in coverage counts");
			require(memberTable==null || !memberTable.next(), "Unexpected extra member coverage row");
		}
		validateSummary(prefix+".summary.tsv", policyTokens, totals, counts, measured);
		if(!bound && baselinePath!=null){writeBaseline();}
		bound=true;
		return counts;
	}

	private void writeBaseline(){
		final ByteBuilder buffer=new ByteBuilder("input_id");
		for(final String name:BASELINE_COLUMNS){buffer.tab().append(name);}
		for(final String name:memberColumns()){buffer.tab().append(name);}
		try(Output output=new Output(baselinePath)){
			output.writer.print(buffer.nl());
			for(int i=0; i<bindings.length; i++){
				BuildNcrnaReference.field(buffer.clear(), bindings[i].id);
				for(int c=0; c<3; c++){buffer.tab().append(baseline[i*3+c]);}
				for(int c=0; c<members.length*2; c++){buffer.tab().append(baselineMembers[i*members.length*2+c]);}
				output.writer.print(buffer.nl());
			}
		}
	}
	private String[] memberColumns(){
		final String[] names=new String[members.length*2];
		for(int m=0; m<members.length; m++){
			names[m*2]="member_"+m+"_TP";
			names[m*2+1]="member_"+m+"_references";
		}
		return names;
	}

	/** Reuses the first round's immutable vectors; this comparison never reruns a caller. */
	private void checkBaseline(final long[] expected, final int offset, final long[] actual){
		assert(offset>=0 && offset<=expected.length-actual.length) : "Each shipped vector occupies its exact flat-array slice";
		if(bound){
			for(int c=0; c<actual.length; c++){
				require(expected[offset+c]==actual[c], "Shipped per-input detection/support vector changed between rounds");
			}
		}else{System.arraycopy(actual, 0, expected, offset, actual.length);}
	}

	private void readMembers(final Table table, final int i, final byte[][] policies,
			final Counts[] counts, final long[] unionTp, final long unionRefs, final long[] root, final long[] row){
		final int input=table.column("input_id"), policy=table.column("policy"), member=table.column("member");
		final int tp=table.column("matched_scored_parents"), refs=table.column("qualified_scored_parents");
		assert(root.length==members.length*2 && row.length==root.length) : "Reusable member scratch follows the bound roster";
		for(int p=0; p<policies.length; p++){
			long sum=0, sumReferences=0;
			for(int m=0; m<members.length; m++){
				require(table.next() && table.matches(inputTokens[i], input) && table.matches(policies[p], policy)
					&& table.matches(memberTokens[m], member), "Missing, duplicate or unexpected input/policy/member coverage row");
				final long matched=count(table, tp), qualified=count(table, refs);
				require(matched<=qualified && matched<=unionTp[p] && qualified<=unionRefs,
					"Member matches/references cannot exceed qualified truth or corresponding union counts");
				sum=Math.addExact(sum, matched);
				sumReferences=Math.addExact(sumReferences, qualified);
				row[m*2]=matched;
				row[m*2+1]=qualified;
				if(p==0){
					root[m*2]=matched;
					root[m*2+1]=qualified;
				}
				require(qualified==root[m*2+1], "Every policy must use the same qualified member-support denominator");
				if(root[m*2]>0 && matched==0){counts[p].memberLastLoss[m]++;}
				if(root[m*2]==0 && matched>0){counts[p].memberFirstGain[m]++;}
			}
			if(p==0){checkBaseline(baselineMembers, i*row.length, row);}
			require(sum>=unionTp[p] && sumReferences>=unionRefs,
				"Every scored union parent must have at least one qualified support member");
		}
	}

	private void validateSummary(final String path, final byte[][] policies, final long[][] totals,
			final Counts[] counts, final Measurement[] measured){
		try(Table table=new Table(path)){
			final int policy=table.column("policy"), lost=table.column("last_loss_inputs"), gained=table.column("first_gain_inputs");
			final int[] columns=columns(table, metrics);
			for(int p=0; p<policies.length; p++){
				require(table.next() && table.matches(policies[p], policy), "Coverage summary must contain each policy exactly once in declared order");
				for(int c=0; c<columns.length; c++){
					require(totals[p][c]==count(table, columns[c]), "Input totals disagree with native policy summary: "+metrics[c]);
				}
				require(counts[p].lastLoss==count(table, lost) && counts[p].firstGain==count(table, gained),
					"Native summary coverage counts disagree with the per-input predicate");
				require(measured[p].tp==totals[p][TP] && measured[p].fp==totals[p][FP]
					&& measured[p].lost==totals[p][LOST] && measured[p].cost==totals[p][COST],
					"Coverage validation must bind the same objective used by selection");
			}
			require(!table.next(), "Unexpected extra coverage summary policy");
		}
	}

	private static int[] columns(final Table table, final String[] names){
		assert(table!=null && names!=null) : "Coverage schema names its required native columns";
		final int[] indices=new int[names.length];
		for(int i=0; i<names.length; i++){indices[i]=table.column(names[i]);}
		return indices;
	}
	static byte[] token(final String name){return BuildNcrnaReference.field(new ByteBuilder(), name).toBytes();}

	/** Saved count fields are unsigned decimal longs. Check before multiplication:
	 * Parse.parseLong's fast path can wrap overflow into a plausible positive count. */
	static long count(final Table table, final int column){
		table.parser.setBounds(column);
		final byte[] line=table.parser.line();
		final int from=table.parser.a(), to=table.parser.b();
		require(from<to, "Empty coverage count in "+table.path);
		long value=0;
		for(int i=from; i<to; i++){
			final int digit=line[i]-'0';
			if(digit<0 || digit>9 || value>(Long.MAX_VALUE-digit)/10){
				throw new IllegalArgumentException("Coverage count must be an unsigned decimal signed-range long in "+table.path);
			}
			value=value*10+digit;
		}
		return value;
	}
	private static void require(final boolean valid, final String message){
		if(!valid){throw new IllegalArgumentException(message);}
	}

	/*--------------------------------------------------------------*/
	/*----------------       Retained Measurements  ----------------*/
	/*--------------------------------------------------------------*/

	/** Immutable after validation. Package writers fill counters before publication;
	 * constructor copies arrays so a test/client cannot mutate a cached point. */
	public static final class Counts{
		public Counts(final long lastLoss_, final long firstGain_, final long[] memberLastLoss_, final long[] memberFirstGain_){
			require(lastLoss_>=0 && firstGain_>=0 && memberLastLoss_!=null && memberFirstGain_!=null
				&& memberLastLoss_.length==memberFirstGain_.length, "Coverage counts need nonnegative union counts and aligned member vectors");
			lastLoss=lastLoss_;
			firstGain=firstGain_;
			memberLastLoss=memberLastLoss_.clone();
			memberFirstGain=memberFirstGain_.clone();
			for(int m=0; m<memberLastLoss.length; m++){
				require(memberLastLoss[m]>=0 && memberFirstGain[m]>=0, "Member coverage counts cannot be negative");
			}
		}
		public long lastLoss(){return lastLoss;}
		public long firstGain(){return firstGain;}
		public long memberLastLoss(final int member){return memberLastLoss[member];}
		public long memberFirstGain(final int member){return memberFirstGain[member];}
		public int members(){return memberLastLoss.length;}
		static boolean same(final Counts a, final Counts b){
			return a==b || a!=null && b!=null && a.lastLoss==b.lastLoss && a.firstGain==b.firstGain
				&& Arrays.equals(a.memberLastLoss, b.memberLastLoss) && Arrays.equals(a.memberFirstGain, b.memberFirstGain);
		}
		long lastLoss, firstGain;
		final long[] memberLastLoss, memberFirstGain;
	}

	private final Binding[] bindings;
	private final String[] metrics;
	private final String[] members;
	private final byte[][] inputTokens, memberTokens;
	private final long[] baseline, baselineMembers;
	private boolean bound;
	private String baselinePath;
	private static final int TP=0, FP=1, FN=2, SCORED=3, MASKED=4, LOST=5, GAINED=6, COST=7;
	private static final String[] METRICS={"TP", "FP", "FN", "scored_references", "masked_references", "lost_baseline_TP", "gained_TP", "kernel_logical_alignments"};
	private static final String[] BASELINE_COLUMNS={"TP", "scored_references", "masked_references"};
}
