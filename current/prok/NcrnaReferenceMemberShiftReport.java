package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import map.IntLongHashMap;
import map.ObjectIntMap;
import map.ObjectMap;
import map.ObjectSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaReferenceGff.Row;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Outputs;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Exact recorded-reference-model joins for source-excluded cobalamin matches.
 * Attribution uses the canonical GFF's recorded reference_model, which may reflect
 * a protected prior. Its reference_evidence is retained; CM scores do not choose
 * or replace that recorded assignment.
 * Adapted from the native endpoint summary/byte-file template. No calling,
 * rematching, coordinate adjustment or inference occurs. All coordinates remain
 * one-based inclusive; signed errors follow transcript direction. Duplicate-check
 * keys cover all selected matches, including excluded sources; joined Match rows
 * are retained only after exclusion. Histogram numbers use primitive maps.
 * @author Brian, Ganyu
 */
public final class NcrnaReferenceMemberShiftReport{

	public static void main(final String[] args){
		final PreParser pre=new PreParser(args, NcrnaReferenceMemberShiftReport.class, false);
		try{
			new NcrnaReferenceMemberShiftReport(pre.args).process();
			pre.outstream.println("NCRNA_REFERENCE_MEMBER_SHIFT_REPORT_COMPLETE");
		}finally{Shared.closeStream(pre.outstream);}
	}

	/** Expected member counts bind this description to an accepted saved grade. */
	public NcrnaReferenceMemberShiftReport(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("bindings")){bindings=b;}
			else if(a.equals("exclude")){exclusions=b;}
			else if(a.equals("expectedall")){expectedAll=counts(b);}
			else if(a.equals("expected")){expected=counts(b);}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown argument: "+arg);}
		}
		input=parser.in1;
		prefix=parser.out1;
		require(input!=null && prefix!=null && bindings!=null && exclusions!=null && expectedAll!=null && expected!=null,
			"Require in=matches bindings=source_map exclude=panel_audit out=prefix expectedall=main,second,third expected=main,second,third");
		require(parser.in2==null && parser.out2==null && !parser.overwrite && !parser.append
			&& parser.maxReads<0 && parser.samplerate==1, "Read complete saved inputs and write fresh reports only");
		final String[] paths={input, bindings, exclusions, prefix+".summary.tsv", prefix+".histogram.tsv", prefix+".loci.tsv", prefix+".population.tsv"};
		require(Tools.testInputFiles(false, true, input, bindings, exclusions)
			&& Tools.testForDuplicateFiles(true, paths)
			&& Tools.testOutputFiles(false, false, false, Arrays.copyOfRange(paths, 3, paths.length)),
			"Inputs and fresh output paths must be distinct and accessible");
		outputs=new FileFormat[4];
		for(int i=0; i<outputs.length; i++){outputs[i]=FileFormat.testOutput(paths[i+3], FileFormat.TXT, null, true, false, false, false);}
		for(int member=0; member<3; member++){
			for(int group=0; group<GROUPS.length; group++){groups[member][group]=new Group();}
		}
	}

	public void process(){
		readBindings();
		readExclusions();
		readMatches();
		for(int member=0; member<3; member++){
			require(selected[member]==expectedAll[member] && retained[member]==expected[member],
				"Saved population differs for "+MEMBERS[member]+": selected="+selected[member]+", retained="+retained[member]);
		}
		joinReferences();
		write();
	}

	/** One input identity names one assembly and one canonical GFF; no taxID joins. */
	private void readBindings(){
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(bindings)){
			RrnaResourceIO.header(in.reader, BINDING_HEADER, bindings);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==5 && unsigned(p, 0)==sources.size(), "Bindings require five fields and dense zero-based source ordinals");
				final String assembly=p.parseString(1), id=p.parseString(2), gff=p.parseString(4);
				require(!assembly.isEmpty() && !assembly.equals(".") && !id.isEmpty() && gff.startsWith("/"),
					"Bindings need an exact assembly, input ID and absolute canonical GFF path");
				for(int i=0; i<id.length(); i++){require(id.charAt(i)>=33 && id.charAt(i)<=126, "Length-prefixed input IDs must be printable ASCII");}
				probe.set(p, 2);
				require(!byInput.contains(probe), "Duplicate input identity in source bindings: "+id);
				final Source source=new Source(id, assembly, gff);
				byInput.put(probe.copy(), source);
				sources.add(source);
			}
		}
		require(!sources.isEmpty(), "The source binding must not be empty");
	}

	/** Mirrors the accepted endpoint report's exact assembly exclusion contract. */
	private void readExclusions(){
		final ObjectSet<String> omitted=new ObjectSet<String>(String.class);
		final ObjectIntMap<String> panels=new ObjectIntMap<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(exclusions)){
			RrnaResourceIO.header(in.reader, EXCLUSION_LABEL, exclusions);
			RrnaResourceIO.header(in.reader, EXCLUSION_HEADER, exclusions);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==5, "Panel audit requires five fields");
				final String panel=p.parseString(0);
				final int previous=panels.get(panel), ordinal=unsigned(p, 1);
				require(panel.startsWith("/") && p.termStartsWith("/", 2)
					&& ordinal==(previous<0 ? 1 : previous+1), "Panel paths and one-based row ordinals must be explicit and dense");
				panels.put(panel, ordinal);
				auditRows++;
				if(p.termEquals("EXCLUDED_ASSEMBLY", 4)){
					final String assembly=p.parseString(3);
					require(!assembly.isEmpty() && !assembly.equals("."), "Excluded rows must name an assembly");
					omitted.add(assembly);
				}else{
					require(p.termEquals("UNMATCHED_MANIFEST_PATH_ALIAS", 4) && p.termEquals(".", 3), "Unknown exclusion status");
					unmatchedAliases++;
				}
			}
		}
		auditPanels=panels.size();
		excludedAssemblies=omitted.size();
		require(auditPanels>0, "An explicit nonempty exclusion audit is required");
		for(final Source source:sources){source.excluded=omitted.contains(source.assembly);}
	}

	/** Validate one-to-one matches before exclusion; never use overlap to infer a CM. */
	private void readMatches(){
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(input)){
			RrnaResourceIO.header(in.reader, MATCH_HEADER, input);
			for(byte[] line=in.reader.nextLine(); line!=null; line=in.reader.nextLine()){
				p.set(line);
				require(p.terms()==9, "Matched rows require the complete nine-column schema");
				int member=-1;
				for(int i=0; i<3; i++){if(p.termEquals(MEMBERS[i], 0)){member=i;}}
				if(member<0){continue;}//Clan rows duplicate member rows; never count them again.
				final int ts=unsigned(p, 3), te=unsigned(p, 4), cs=unsigned(p, 5), ce=unsigned(p, 6);
				final boolean minus=p.termEquals("-", 2);
				require((minus || p.termEquals("+", 2)) && ts>=1 && te>=ts && cs>=1 && ce>=cs,
					"Matches require explicit strand and nonempty one-based inclusive coordinates");
				final long overlap=Math.max(0L, (long)Math.min(te, ce)-Math.max(ts, cs)+1);
				require(2*overlap>=Math.max((long)te-ts+1, (long)ce-cs+1), "Saved match violates reciprocal-half overlap");
				final int d5=minus ? te-ce : cs-ts, d3=minus ? ts-cs : ce-te;
				require(signed(p, 7)==d5 && signed(p, 8)==d3, "Saved signed errors disagree with transcript-oriented coordinates");
				p.setBounds(1);
				final int from=p.a(), limit=p.b();
				int colon=from;
				while(colon<limit && line[colon]>='0' && line[colon]<='9'){colon++;}
				require(colon>from && colon<limit && line[colon]==':', "Matched seqid must use length:input_id/contig");
				final int length=NcrnaNumbers.unsignedInt(line, from, colon);
				final long separator=(long)colon+1+length;
				require(length>0 && separator<limit-1 && line[(int)separator]=='/', "Source-ID length must end at the contig separator");
				probe.set(line, colon+1, (int)separator);
				final Source source=byInput.get(probe);
				require(source!=null, "Matched input ID is absent from the source binding");
				final String contig=new String(line, (int)separator+1, limit-(int)separator-1, StandardCharsets.UTF_8);
				key(contig, minus ? '-' : '+', ts, te);
				final Target truth=probe.copy();
				require(source.truth.add(truth), "A clan reference was consumed by multiple member matches");
				key(contig, minus ? '-' : '+', cs, ce);
				require(source.calls.add(probe.copy()), "A clan call was consumed by multiple member matches");
				selected[member]++;
				if(source.excluded){continue;}
				retained[member]++;
				final Match match=new Match(source, member, contig, minus ? '-' : '+', ts, te, cs, ce, d5, d3);
				source.wanted.put(truth, match);
				matches.add(match);
			}
		}
		assert(matches.size()==retained[0]+retained[1]+retained[2]) : "Only retained member rows enter the reference join";
	}

	/** Canonical codec validates file/input framing, attributes, escapes and locus IDs. */
	private void joinReferences(){
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		for(final Source source:sources){
			if(source.wanted.isEmpty()){continue;}
			referenceFiles++;
			for(final Row row:codec.read(source.gff, source.id)){
				key(row.locus.seqid, row.locus.strand, row.locus.start, row.locus.stop);
				final Match match=source.wanted.get(probe);
				if(match==null){continue;}
				match.exactRows++;
				if(row.mask || !row.locus.family.equals("cobalamin")){continue;}
				match.scoredClanRows++;
				// Never select the last of multiple exact rows as the winning CM.
				match.referenceModel=match.scoredClanRows==1 ? row.model : null;
				match.referenceId=match.scoredClanRows==1 ? row.id : null;
				match.referenceEvidence=match.scoredClanRows==1 ? row.evidence : match.referenceEvidence+"\n"+row.evidence;
			}
		}
		for(final Match match:matches){
			int group;
			if(match.exactRows==0){group=3;}
			else if(match.scoredClanRows==0){group=7;}
			else if(match.scoredClanRows>1){group=6;}
			else{
				group=4;
				if(match.referenceModel!=null){
					if(match.referenceModel.indexOf(',')>=0){group=5;}
					for(int i=0; i<3; i++){if(MEMBERS[i].equals(match.referenceModel)){group=i;}}
				}
			}
			match.group=group;
			groups[match.member][group].add(match.d5, match.d3);
			groups[match.member][8].add(match.d5, match.d3);
		}
		for(int member=0; member<3; member++){
			long sum=0;
			for(int group=0; group<8; group++){sum+=groups[member][group].n;}
			require(sum==retained[member] && groups[member][8].n==retained[member], "Resolved, ambiguous and unresolved groups must conserve each member population");
		}
	}

	/** Reuses a byte probe; only retained match keys copy these bytes. */
	private void key(final String contig, final char strand, final int start, final int stop){
		assert(contig!=null && !contig.isEmpty() && start>=1 && stop>=start) : "Exact keys require a literal nonempty contig and inclusive interval";
		buffer.clear().append(contig).tab().append(strand).tab().append(start).tab().append(stop);
		probe.set(buffer.array, 0, buffer.length());
	}

	private void write(){
		try(Outputs owners=new Outputs(4)){
			final ByteStreamWriter summary=owners.open(outputs[0]);
			final ByteStreamWriter histogram=owners.open(outputs[1]);
			final ByteStreamWriter loci=owners.open(outputs[2]);
			final ByteStreamWriter population=owners.open(outputs[3]);
			summary.print("# Empty groups have no distribution and are omitted; population.tsv retains selected/retained/excluded member counts.\n");
			summary.print("query_member\trecorded_reference_group\tend\tmatched\twithin10\twithin10_pct\tminimum\tp10\tp25\tp50\tp75\tp90\tmaximum\tcentral80_width_bp\tiqr_bp\tlowest_mode_bp\tmode_loci\tmode_pct\ttied_modal_bins\n");
			histogram.print("query_member\trecorded_reference_group\tend\tsigned_offset_bp\tloci\n");
			population.print("query_member\tselected\tretained\texcluded\tbinding_sources\texcluded_assemblies\taudit_panels\taudit_rows\tunmatched_alias_rows\treference_files\n");
			for(int member=0; member<3; member++){
				population.print(buffer.clear().append(MEMBERS[member]).tab().append(selected[member]).tab().append(retained[member])
					.tab().append(selected[member]-retained[member]).tab().append(sources.size()).tab().append(excludedAssemblies)
					.tab().append(auditPanels).tab().append(auditRows).tab().append(unmatchedAliases).tab().append(referenceFiles).nl());
				for(int group=0; group<GROUPS.length; group++){
					for(int end=0; end<2; end++){distribution(summary, histogram, member, group, end);}
				}
			}
			loci.print("query_member\tsource\tassembly\tcontig\tstrand\ttruth_start\ttruth_stop\tcall_start\tcall_stop\td5\td3\trecorded_reference_group\texact_gff_rows\tscored_clan_rows\trecorded_reference_model\treference_evidence\treference_id\treference_gff\n");
			for(final Match match:matches){
				buffer.clear().append(MEMBERS[match.member]).tab().append(match.source.id).tab().append(match.source.assembly).tab();
				BuildNcrnaReference.field(buffer, match.contig).tab().append(match.strand).tab().append(match.ts).tab().append(match.te)
					.tab().append(match.cs).tab().append(match.ce).tab().append(match.d5).tab().append(match.d3)
					.tab().append(GROUPS[match.group]).tab().append(match.exactRows).tab().append(match.scoredClanRows).tab();
				BuildNcrnaReference.field(buffer, match.referenceModel==null ? "." : match.referenceModel).tab();
				BuildNcrnaReference.field(buffer, match.referenceEvidence==null ? "." : match.referenceEvidence).tab();
				BuildNcrnaReference.field(buffer, match.referenceId==null ? "." : match.referenceId).tab();
				BuildNcrnaReference.field(buffer, match.source.gff).nl();
				loci.print(buffer);
			}
		}
	}

	/** Nearest-rank quantiles and central spread are descriptive, not a correction. */
	private void distribution(final ByteStreamWriter summary, final ByteStreamWriter histogram, final int member, final int group, final int end){
		final Group stats=groups[member][group];
		if(stats.n==0){return;}// No distribution exists; omit the row rather than invent its statistics.
		final IntLongHashMap values=stats.signed[end];
		final int[] keys=values.toArray();
		Arrays.sort(keys);
		final String label=end==0 ? "5prime" : "3prime";
		long total=0, within=0, modeN=0;
		int mode=0, tied=0;
		for(final int key:keys){
			final long count=values.get(key);
			total+=count;
			if(key>=-10 && key<=10){within+=count;}
			if(count>modeN){
				modeN=count;
				mode=key;
				tied=1;
			}else if(count==modeN){tied++;}
			histogram.print(buffer.clear().append(MEMBERS[member]).tab().append(GROUPS[group]).tab().append(label).tab().append(key).tab().append(count).nl());
		}
		require(total==stats.n, "A signed histogram must account for every joined locus once");
		buffer.clear().append(MEMBERS[member]).tab().append(GROUPS[group]).tab().append(label).tab().append(total).tab().append(within);
		{
			final int[] quantiles=new int[5];
			long cumulative=0;
			int index=0;
			for(int q=0; q<5; q++){
				final int percent=PERCENTILES[q];
				final long rank=total/100*percent+(total%100*percent+99)/100;
				while(cumulative<rank){cumulative+=values.get(keys[index++]);}
				quantiles[q]=keys[index-1];
			}
			buffer.tab().append(100.0*within/total, 6).tab().append(keys[0]);
			for(final int quantile:quantiles){buffer.tab().append(quantile);}
			buffer.tab().append(keys[keys.length-1]).tab().append((long)quantiles[4]-quantiles[0])
				.tab().append((long)quantiles[3]-quantiles[1]).tab().append(mode).tab().append(modeN)
				.tab().append(100.0*modeN/total, 6).tab().append(tied);
		}
		summary.print(buffer.nl());
	}

	private static long[] counts(final String text){
		require(text!=null, "Expected counts must explicitly name main, second and third members");
		final String[] fields=text.split(",");
		require(fields.length==3, "Expected counts use RF00174,RF01482,RF01689 order");
		final long[] values=new long[3];
		for(int i=0; i<3; i++){
			values[i]=Long.parseLong(fields[i]);
			require(values[i]>=0, "Expected populations cannot be negative");
		}
		return values;
	}
	private static int unsigned(final LineParser1 p, final int field){
		p.setBounds(field);
		return NcrnaNumbers.unsignedInt(p.line(), p.a(), p.b());
	}
	private static int signed(final LineParser1 p, final int field){
		p.setBounds(field);
		final boolean negative=p.a()<p.b() && p.line()[p.a()]=='-';
		final int value=NcrnaNumbers.unsignedInt(p.line(), p.a()+(negative ? 1 : 0), p.b());
		return negative ? -value : value;
	}
	private static void require(final boolean valid, final String reason){if(!valid){throw new IllegalArgumentException(reason);}}

	private static final class Source{
		Source(final String id_, final String assembly_, final String gff_){
			id=id_;
			assembly=assembly_;
			gff=gff_;
		}
		final String id, assembly, gff;
		boolean excluded;
		final ObjectSet<Target> truth=new ObjectSet<Target>(Target.class), calls=new ObjectSet<Target>(Target.class);
		final ObjectMap<Target, Match> wanted=new ObjectMap<Target, Match>(Target.class, Match.class);
	}
	private static final class Match{
		Match(final Source source_, final int member_, final String contig_, final char strand_, final int ts_, final int te_,
				final int cs_, final int ce_, final int d5_, final int d3_){
			source=source_;
			member=member_;
			contig=contig_;
			strand=strand_;
			ts=ts_;
			te=te_;
			cs=cs_;
			ce=ce_;
			d5=d5_;
			d3=d3_;
		}
		final Source source;
		final String contig;
		final int member, ts, te, cs, ce, d5, d3;
		final char strand;
		int exactRows, scoredClanRows, group;
		String referenceModel, referenceId, referenceEvidence;
	}
	private static final class Group{
		void add(final int five, final int three){
			assert(n<Long.MAX_VALUE) : "Each joined locus increments each endpoint histogram once";
			n++;
			signed[0].increment(five, 1);
			signed[1].increment(three, 1);
		}
		long n;
		final IntLongHashMap[] signed={new IntLongHashMap(16), new IntLongHashMap(16)};
	}

	private String input, bindings, exclusions, prefix;
	private long[] expectedAll, expected;
	private FileFormat[] outputs;
	private final long[] selected=new long[3], retained=new long[3];
	private final Group[][] groups=new Group[3][GROUPS.length];
	private final ArrayList<Source> sources=new ArrayList<Source>();
	private final ArrayList<Match> matches=new ArrayList<Match>();
	private final ObjectMap<Target, Source> byInput=new ObjectMap<Target, Source>(Target.class, Source.class);
	private final Target probe=new Target();
	private final ByteBuilder buffer=new ByteBuilder();
	private int auditPanels, excludedAssemblies, referenceFiles;
	private long auditRows, unmatchedAliases;
	private static final String[] MEMBERS={"RF00174", "RF01482", "RF01689"};
	private static final String[] GROUPS={"RF00174", "RF01482", "RF01689", "UNRESOLVED_NO_EXACT_REFERENCE",
		"UNRESOLVED_MODEL", "AMBIGUOUS_MODEL", "AMBIGUOUS_EXACT_REFERENCE", "UNRESOLVED_NOT_SCORED", "ALL_REFERENCE_MODELS"};
	private static final int[] PERCENTILES={10,25,50,75,90};
	private static final String BINDING_HEADER="input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff";
	private static final String EXCLUSION_LABEL="# HELD OUT (panel genomes excluded from tables and training): exclusion bindings; not a fitted result";
	private static final String EXCLUSION_HEADER="panel\trow1\tfasta\tassembly\tstatus";
	private static final String MATCH_HEADER="#family\tseqid\tstrand\ttruth_start\ttruth_stop\tcall_start\tcall_stop\td5\td3";
}
