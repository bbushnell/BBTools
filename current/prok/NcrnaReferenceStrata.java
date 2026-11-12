package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import map.LongHashSet;
import map.LongIntMap;
import map.ObjectMap;
import parse.LineParser1;
import prok.NcrnaReferenceEvidence.Model;
import prok.NcrnaReferenceEvidence.Roster;
import prok.NcrnaReferenceGeometry.Hit;
import prok.NcrnaReferenceGeometry.InputState;
import prok.NcrnaReferenceGeometry.Pair;
import prok.NcrnaReferenceGeometry.Target;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import structures.ByteBuilder;
import structures.IntList;
import structures.LongList;

/** Descriptive event/source strata over saved reference evidence.
 * Reuses the native geometry loader; opens no genome, GFF or CM payloads.
 * Flags are nonexclusive observations, never biological admission or exclusion.
 * RFAM claims are source assertions, not independent model confirmations.
 * Event attribute flags apply only to directly overlapping RNA-type source rows;
 * context/component/other flags remain in qualified links and source rows. No
 * attribute is inherited through Parent links. Ancestors remain separate rows.
 * One calling thread owns all state. Parsers and temporary byte buffers are reused;
 * only source keys, compact row state and native model references are retained.
 * @author Brian Bushnell, Raiden
 */
final class NcrnaReferenceStrata{

	/*--------------------------------------------------------------*/
	/*----------------         Initialization       ----------------*/
	/*--------------------------------------------------------------*/

	static void process(final String input, final String attributes, final String attributeEvidence, final String output,
		final NcrnaEvidenceCatalog catalog, final Roster roster){
		final NcrnaReferenceGeometry geometry=NcrnaReferenceGeometry.load(input, catalog, roster, true);
		load(input, attributes, attributeEvidence, geometry, catalog, roster).write(output);
	}

	/** Reuses the saved source join without materializing and rereading a full report. */
	static NcrnaReferenceStrata load(final String input, final String attributes, final String attributeEvidence,
		final NcrnaReferenceGeometry geometry, final NcrnaEvidenceCatalog catalog, final Roster roster){
		final NcrnaReferenceStrata strata=new NcrnaReferenceStrata(geometry, roster, catalog.size());
		strata.readHistory(input+".historical.tsv.gz");
		strata.readRelations(input+".relations.tsv.gz");
		final NcrnaReferenceProfiles captured=attributeEvidence==null ? geometry.profiles()
			: NcrnaReferenceProfiles.sourceClaims(attributeEvidence, catalog);
		if(attributeEvidence!=null){
			require(captured.observed()==catalog.size() && strata.profiles.observed()==catalog.size(),
				"Cross-generation attribute reuse requires the same complete supplied catalog");
			for(int i=0; i<catalog.size(); i++){
				final String id=catalog.get(i).id;
				require(captured.sameSources(strata.profiles, id),
					"Captured source claims/states/accounting differ from the new source join for input "+id);
			}
		}
		strata.readAttributes(attributes, captured.evidenceState());
		strata.readOverlaps(input+".overlaps.tsv.gz");
		strata.addGeometry();
		return strata;
	}
	/** Borrowed immutable-after-load data for the single-threaded reference constructor. */
	int[] flags(){return flags;}
	ArrayList<SourceRow> sourceRows(){return rows;}
	IntList linkedSources(){return linkSources;}
	IntList linkedHits(){return linkHits;}

	private NcrnaReferenceStrata(final NcrnaReferenceGeometry geometry_, final Roster roster_, final int catalogInputs_){
		geometry=geometry_;
		roster=roster_;
		catalogInputs=catalogInputs_;
		profiles=geometry.profiles();
		flags=new int[geometry.hits().length];
		sameRows=new int[flags.length];
		oppositeRows=new int[flags.length];
		originCounts=new int[4][flags.length];
		rfamSeen=new boolean[flags.length];
	}

	/*--------------------------------------------------------------*/
	/*----------------         Saved sources        ----------------*/
	/*--------------------------------------------------------------*/

	private void readHistory(final String path){
		final LineParser1 p=new LineParser1('\t');
		final ObjectMap<Target,Target> targets=new ObjectMap<Target,Target>(Target.class, Target.class);
		final Target probe=new Target();
		SourceFile source=null;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, profiles.batchEvidence() ? BuildNcrnaReference.BATCH_HISTORICAL_HEADER : BuildNcrnaReference.HISTORICAL_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==(profiles.batchEvidence() ? 20 : 19), "Historical strata schema differs from its evidence generation");
				source=source(p, 0, 1, source, true);
				final long line=p.parseLong(2);
				require(line>0 && !source.lines.contains(line), "Repeated or invalid input/GFF/line in historical strata input");
				probe.set(p, 3);
				Target target=targets.get(probe);
				if(target==null){
					target=probe.copy();
					targets.put(target, target);
				}
				probe.set(p, 4);
				Target type=targets.get(probe);
				if(type==null){
					type=probe.copy();
					targets.put(type, type);
				}
				final SourceRow item=new SourceRow(source, line, target, type, rows.size(), p);
				source.lines.put(line, rows.size());
				rows.add(item);
			}
		}
		require(rows.size()==profiles.retained(), "Source strata must retain every saved historical row");
		for(final SourceInput input:inputs.values()){
			if(input==null){continue;}
			long count=0;
			for(final SourceFile file:input.files.values()){if(file!=null){count+=file.lines.size();}}
			require(count==profiles.retained(input.state.genome.id), "Historical row count differs for input "+input.state.genome.id);
		}
	}

	private void readRelations(final String path){
		final LineParser1 p=new LineParser1('\t');
		SourceFile source=null;
		String lastAccession=null;
		Model lastModel=null;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, BuildNcrnaReference.RELATION_HEADER, path);
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==7, "Saved source relations need seven fields");
				if(p.termEquals("CYCLE", 3)){continue;}
				source=source(p, 0, 1, source, false);
				final SourceRow item=row(source, p.parseLong(2));
				if(p.termEquals("RFAM", 3)){
					require(item.modelsSeen<item.models.length, "Too many RFAM claims for historical row");
					if(lastAccession==null || !p.termEquals(lastAccession, 5)){
						lastAccession=text(p, 5);
						lastModel=roster.model(lastAccession);
					}
					item.models[item.modelsSeen++]=lastModel;
					if(lastModel==null){item.flags|=UNKNOWN_RFAM;}
				}else{
					require(p.termEquals("PARENT", 3), "Unknown saved source relationship kind");
					item.parentsSeen++;
					if(p.termEquals("DANGLING", 6)){item.flags|=DANGLING_PARENT;}
				}
			}
		}
		for(final SourceRow item:rows){
			if(item.modelsSeen!=item.models.length || item.parentsSeen!=item.parents){
				throw new IllegalArgumentException("Saved source relationship counts do not match history at "
					+new String(item.source.path, StandardCharsets.UTF_8)+":"+item.line);
			}
		}
	}

	/** Reuses captured claims under their original completed evidence status.
	 * Explicit cross-generation reuse validates the catalog/source population first;
	 * this join then requires every exact input/GFF/line key once. Accepted producer
	 * provenance identifies the frozen source corpus; no live GFF reread or new claim
	 * about uncaptured historical bytes is implied by a matching row count. */
	private void readAttributes(final String prefix, final String capturedState){
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(prefix+".attribute_summary.tsv")){
			RrnaResourceIO.header(in.reader, ATTRIBUTE_SUMMARY_HEADER, prefix);
			final byte[] row=in.reader.nextLine();
			require(row!=null, "Attribute capture must have a completed summary");
			p.set(row);
			require(p.terms()==9 && p.parseInt(0)==profiles.observed() && p.parseInt(1)==catalogInputs && p.parseLong(2)==rows.size()
				&& p.parseLong(3)==rows.size() && p.termEquals(capturedState, 7)
				&& p.termEquals("CURRENT_SOURCE_ATTRIBUTES_NOT_ADJUDICATED", 8) && in.reader.nextLine()==null,
				"Attribute capture population/status differs from the saved evidence");
		}
		SourceFile source=null;
		final Claims claims=new Claims();
		long count=0;
		try(Input in=new Input(prefix+".attributes.tsv.gz")){
			RrnaResourceIO.header(in.reader, "input_id\tgff\tline\traw_attributes", prefix);
			for(byte[] raw=in.reader.nextLine(); raw!=null; raw=in.reader.nextLine()){
				p.set(raw);
				require(p.terms()==4, "Captured attributes require four fields");
				source=source(p, 0, 1, source, false);
				final SourceRow item=row(source, p.parseLong(2));
				require(!item.attributesSeen, "Repeated captured input/GFF/line key");
				item.attributesSeen=true;
				p.setBounds(3);
				item.flags|=claims.read(raw, p.a(), p.b());
				count++;
			}
		}
		require(count==rows.size(), "Captured attributes must cover every historical key exactly once");
	}

	private SourceFile source(final LineParser1 p, final int inputField, final int fileField,
		final SourceFile cached, final boolean create){
		if(cached!=null && p.termEquals(cached.input.state.genome.id, inputField) && same(p, fileField, cached.path)){return cached;}
		final String id=text(p, inputField);
		SourceInput input=inputs.get(id);
		if(input==null && create){
			final InputState state=geometry.input(id);
			require(state!=null && state.selected, "Source key must belong to the saved observed population: "+id);
			input=new SourceInput(state);
			inputs.put(id, input);
		}
		require(input!=null, "Unknown input in saved source join: "+id);
		final String encoded=text(p, fileField);
		SourceFile source=input.files.get(encoded);
		if(source==null && create){
			final ByteBuilder buffer=new ByteBuilder();
			boolean claimed=false;
			for(int i=0; i<input.state.genome.memberships(); i++){
				BuildNcrnaReference.field(buffer.clear(), input.state.genome.membership(i).gff);
				claimed|=same(p, fileField, buffer.toBytes());
			}
			require(claimed, "Saved GFF must be claimed by this exact catalog input");
			source=new SourceFile(input, p.parseByteArray(fileField));
			input.files.put(encoded, source);
		}
		require(source!=null, "Unknown GFF in saved source join: "+encoded);
		return source;
	}

	private SourceRow row(final SourceFile source, final long line){
		final int index=source.lines.get(line);
		require(index>=0, "Saved relationship/attribute/overlap lacks its historical source key");
		return rows.get(index);
	}

	/*--------------------------------------------------------------*/
	/*----------------          Event links         ----------------*/
	/*--------------------------------------------------------------*/

	private void readOverlaps(final String path){
		final LineParser1 p=new LineParser1('\t');
		final LongHashSet seen=new LongHashSet(16, 0.7f);
		SourceFile source=null;
		long count=0;
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, BuildNcrnaReference.OVERLAP_HEADER, path);
			for(byte[] raw=in.reader.nextLine(); raw!=null; raw=in.reader.nextLine()){
				p.set(raw);
				require(p.terms()==9, "Saved overlap needs nine fields");
				if(p.termEquals("MODEL_MODEL", 1)){continue;}
				require(p.termEquals("HISTORICAL_MODEL", 1), "Unknown overlap kind");
				source=source(p, 0, 2, source, false);
				final SourceRow item=row(source, p.parseLong(3));
				final InputState input=source.input.state;
				final int index=input.find(p, 4, p.parseLong(5));
				require(index>=0, "Historical overlap names no retained native event");
				require(seen.add(((long)item.index<<32)|(index&0xffffffffL)), "Repeated historical/native overlap pair");
				final Hit hit=geometry.hits()[index];
				final long overlap=Math.max(0L, Math.min(item.stop, hit.stop)-Math.max(item.start, hit.start)+1);
				require(item.assessable && item.target.equals(hit.target) && overlap>0 && overlap==p.parseLong(6),
					"Historical link geometry disagrees with its saved source/native rows");
				if(!hit.valid()){
					require(p.termEquals("UNASSESSABLE_HIT_RANGE", 7), "Invalid hit cannot supply source support");
					item.invalidSeen++;
				}else if(hit.strand!=item.strand){
					require(p.termEquals("OPPOSITE_STRAND_GEOMETRY", 7), "Opposite-strand source overlap must remain separate");
					item.oppositeSeen++;
					oppositeRows[index]++;
					flags[index]|=OPPOSITE_SOURCE;
				}else{
					require(p.termEquals("SAME_STRAND_CANDIDATE", 7), "Same-strand overlap must retain its source relation");
					item.sameSeen++;
					sameRows[index]++;
					originCounts[item.kind][index]++;
					// A coding/context row's pseudo flag is not an RNA pseudo claim.
					if(item.kind==RNA){flags[index]|=item.flags;}
					else{flags[index]|=item.flags&UNKNOWN_RFAM;}
					rfamSeen[index]|=item.models.length>0;
					flags[index]|=compareClaims(item, hit)&~NO_SOURCE_RFAM;
				}
				linkSources.add(item.index);
				linkHits.add(index);
				count++;
			}
		}
		require(count==profiles.historicalPairs(), "Historical event links must conserve the saved denominator");
		for(final SourceRow item:rows){
			if(item.sameSeen!=item.sameExpected || item.oppositeSeen!=item.oppositeExpected || item.invalidSeen!=item.invalidExpected){
				throw new IllegalArgumentException("Historical link counts disagree at "
					+new String(item.source.path, StandardCharsets.UTF_8)+":"+item.line);
			}
		}
	}

	/** A literal RFAM comparison is retained even for regulatory or tagged gene rows.
	 * Its spatial/strand association is a separate field, not an admission assertion.
	 */
	private int compareClaims(final SourceRow source, final Hit hit){
		if(source.models.length==0){return NO_SOURCE_RFAM;}
		final Model own=roster.query(hit.query);
		int result=0;
		for(final Model model:source.models){
			if(model==null){result|=UNKNOWN_RFAM;}
			else{
				if(model.relatedClass()){result|=RELATED_SOURCE_CLASS;}
				result|=model==own ? SAME_ACCESSION_SOURCE : model.scope.family.equals(own.scope.family) ? SAME_FAMILY_SOURCE : OTHER_FAMILY_SOURCE;
			}
		}
		return result;
	}

	private void addGeometry(){
		for(final Hit hit:geometry.hits()){
			int bits=hit.to-hit.from+1==roster.query(hit.query).length ? FULL : SHORT;
			bits|=hit.truncation==0 ? UNTRUNCATED : TRUNCATED;
			if(!hit.valid()){bits|=INVALID;}
			if(hit.left || hit.right){bits|=EDGE;}
			if(hit.valid()){
				if(sameRows[hit.index]==0){bits|=NO_SAME_SOURCE;}
				else if(!rfamSeen[hit.index]){bits|=NO_SOURCE_RFAM;}
			}
			flags[hit.index]|=bits;
		}
		for(final Pair pair:geometry.pairs()){
			final Hit a=pair.a, b=pair.b;
			if(!a.valid() || !b.valid()){continue;}
			int bits=0;
			if(a.strand!=b.strand){bits=OPPOSITE_MODEL;}
			else if(a.query!=b.query){
				bits=OTHER_MODEL;
				if(!roster.query(a.query).scope.family.equals(roster.query(b.query).scope.family)){bits|=OTHER_FAMILY_MODEL;}
				else{bits|=SAME_FAMILY_MODEL;}
			}else{bits=SAME_MODEL;}
			flags[a.index]|=bits;
			flags[b.index]|=bits;
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Output            ----------------*/
	/*--------------------------------------------------------------*/

	private void write(final String prefix){
		final ByteBuilder b=new ByteBuilder();
		final int n=roster.queried(), columns=NAMES.length+1;
		final long[][] events=new long[n+1][columns], groups=new long[n+1][columns];
		final int[] groupFlags=new int[n+1];
		final boolean[] groupSeen=new boolean[n+1];
		final JointCounts joint=new JointCounts();
		int current=0;
		try(Output out=new Output(prefix+".event_strata.tsv.gz")){
			out.writer.println("input_id\tnative_table\tnative_line\tgroup_id\taccession\tmodel_span\tmodel_length\tnative_truncation\tleft_contig_edge\tright_contig_edge\tsame_strand_source_rows\topposite_strand_source_rows\trna_rows\tcomponent_rows\tcontext_rows\tother_rows\tbits\tGA_bits\tGA_margin_bits\tflags\tstatus");
			for(final Hit hit:geometry.ordered()){
				if(hit.group!=current){
					if(current>0){finishGroup(groups, groupFlags, groupSeen);}
					current=hit.group;
				}
				final int bits=flags[hit.index];
				count(events[hit.query], bits);
				count(events[n], bits);
				joint.add(hit, hit.query);
				joint.add(hit, n);
				if(hit.group>0){
					groupSeen[hit.query]=groupSeen[n]=true;
					groupFlags[hit.query]|=bits;
					groupFlags[n]|=bits;
				}
				b.clear().append(hit.input.genome.id).tab().append(hit.table).tab().append(hit.line).tab();
				if(hit.group==0){b.append("NA");}else{b.append(hit.group);}
				b.tab().append(roster.query(hit.query).accession).tab().append(hit.to-hit.from+1).tab().append(roster.query(hit.query).length).tab();
				b.append(hit.truncation).tab().append(hit.left).tab().append(hit.right).tab().append(sameRows[hit.index]).tab().append(oppositeRows[hit.index]).tab();
				for(final int[] counts:originCounts){b.append(counts[hit.index]).tab();}
				b.append(hit.bits, 6).tab().append(roster.query(hit.query).ga, 6).tab().append(hit.bits-roster.query(hit.query).ga, 6).tab();
				out.writer.print(names(b, bits).append("\tEVIDENCE_NOT_ADJUDICATED\n"));
			}
			if(current>0){finishGroup(groups, groupFlags, groupSeen);}
		}
		require(events[n][0]==geometry.hits().length && groups[n][0]==geometry.groups(), "Strata must conserve every native event and review group");
		try(Output out=new Output(prefix+".source_strata.tsv.gz")){
			out.writer.println("input_id\tgff\tline\tsource_type\tstart\tstop\tstrand\tsource_kind\tselection_bits\tsame_strand_hits\topposite_strand_hits\tunassessable_hits\tflags\tstatus");
			for(final SourceRow item:rows){
				b.clear().append(item.source.input.state.genome.id).tab().append(item.source.path).tab().append(item.line).tab();
				b.append(item.type.bytes).tab().append(item.start).tab().append(item.stop).tab().append(item.strand).tab();
				b.append(KINDS[item.kind]).tab().append(item.selection).tab();
				b.append(item.sameSeen).tab().append(item.oppositeSeen).tab().append(item.invalidSeen).tab();
				out.writer.print(names(b, item.flags).append("\tSOURCE_CLAIMS_NOT_ADJUDICATED\n"));
			}
		}
		try(Output out=new Output(prefix+".source_links.tsv.gz")){
			out.writer.println("input_id\tgff\tline\tnative_table\tnative_line\tsource_type\tsource_kind\tassociation\trfam_comparison\tsource_flags");
			for(int i=0; i<linkSources.size; i++){
				final SourceRow item=rows.get(linkSources.get(i));
				final Hit hit=geometry.hits()[linkHits.get(i)];
				b.clear().append(hit.input.genome.id).tab().append(item.source.path).tab().append(item.line).tab();
				b.append(hit.table).tab().append(hit.line).tab().append(item.type.bytes).tab().append(KINDS[item.kind]).tab();
				b.append(!hit.valid() ? "UNASSESSABLE" : hit.strand==item.strand ? "SAME_STRAND_OVERLAP" : "OPPOSITE_STRAND_OVERLAP").tab();
				names(b, compareClaims(item, hit)).tab();
				out.writer.print(names(b, item.flags).nl());
			}
		}
		try(Output out=new Output(prefix+".strata_counts.tsv")){
			out.writer.println("accession\tstratum\tevents\treview_groups\tstatus");
			for(int model=0; model<=n; model++){
				for(int flag=0; flag<columns; flag++){
					b.clear().append(model==n ? "ALL_QUERIED" : roster.query(model).accession).tab().append(flag==0 ? "ALL" : NAMES[flag-1]);
					out.writer.print(b.tab().append(events[model][flag]).tab().append(groups[model][flag]).append("\tNONEXCLUSIVE_NONADDITIVE_EVIDENCE_NOT_GENES\n"));
				}
			}
		}
		joint.write(prefix+".joint_strata.tsv");
		try(Output out=new Output(prefix+".strata_summary.tsv")){
			out.writer.println("inputs_observed\tcatalog_inputs\traw_hit_events\treview_groups\tretained_source_rows\thistorical_model_pairs\tevidence_status\tstatus");
			b.clear().append(profiles.observed()).tab().append(catalogInputs).tab().append(events[n][0]).tab().append(groups[n][0]).tab().append(rows.size()).tab();
			out.writer.print(b.append(profiles.historicalPairs()).tab().append(profiles.evidenceState()).append("\tSTRATA_COMPLETE_NOT_ADJUDICATED\n"));
		}
	}

	/** Exact event membership and distinct group membership in each displayed bin.
	 * Groups can occupy several bins/accessions: their counts are not additive.
	 * The numeric bins are display intervals, never admission thresholds.
	 */
	private final class JointCounts{
		void add(final Hit hit, final int query){
			if(hit.group!=currentGroup){
				assert(hit.group>currentGroup) : "Geometry.ordered supplies nondecreasing review groups for one-group scratch reuse";
				finishGroup();
				currentGroup=hit.group;
			}
			final Model model=roster.query(hit.query);
			final int span=hit.to-hit.from+1, bits=flags[hit.index];
			final int coverage=coverageBin(span, model.length);
			final double difference=hit.bits-model.ga;
			require(Double.isFinite(difference), "Displayed GA-margin intervals require a finite score difference");
			final int margin=marginBin(difference);
			final int source;
			if(!hit.valid()){source=0;}
			else if(sameRows[hit.index]==0){source=profiles.availableClaims(hit.input.genome.id)==0 ? 1 : 2;}
			else{source=originCounts[RNA][hit.index]>0 ? 3 : 4;}
			final int claims=((bits&SAME_ACCESSION_SOURCE)!=0 ? 1 : 0)|((bits&SAME_FAMILY_SOURCE)!=0 ? 2 : 0)
				|((bits&OTHER_FAMILY_SOURCE)!=0 ? 4 : 0)|((bits&UNKNOWN_RFAM)!=0 ? 8 : 0)|((bits&RELATED_SOURCE_CLASS)!=0 ? 16 : 0);
			final int alternatives=((bits&SAME_MODEL)!=0 ? 1 : 0)|((bits&SAME_FAMILY_MODEL)!=0 ? 2 : 0)
				|((bits&OTHER_FAMILY_MODEL)!=0 ? 4 : 0)|((bits&OPPOSITE_MODEL)!=0 ? 8 : 0);
			long key=query;
			key=key*COVERAGE_LABELS.length+coverage;
			key=key*4+hit.truncation;
			key=key*MARGIN_LABELS.length+margin;
			key=key*SOURCE_LABELS.length+source;
			key=key*32+claims;
			// Source RNA claim bits9..16 are contiguous; preserve every contradiction.
			key=key*256+((bits>>>9)&255);
			key=key*16+alternatives;
			key=key*2+(hit.left || hit.right ? 1 : 0);
			events.increment(key);
			if(hit.group>0){groupKeys.add(key);}
		}
		private void finishGroup(){
			Arrays.sort(groupKeys.array, 0, groupKeys.size);
			for(int i=0; i<groupKeys.size; i++){
				if(i==0 || groupKeys.array[i]!=groupKeys.array[i-1]){groups.increment(groupKeys.array[i]);}
			}
			groupKeys.clear();
		}
		void write(final String path){
			finishGroup();
			final long[] keys=events.toArray();
			Arrays.sort(keys);
			final ByteBuilder b=new ByteBuilder();
			try(Output out=new Output(path)){
				out.writer.println("accession\tmodel_coverage\tnative_truncation\tGA_margin_bits\tsource_association\trfam_claim_comparison\trna_claim_flags\tsame_model_overlap\tsame_family_CM_alternative\tother_family_CM_alternative\topposite_strand_CM_overlap\tphysical_contig_edge\tevents\treview_groups\tstatus");
				for(final long key:keys){
					long rest=key;
					final int edge=(int)(rest%2);
					rest/=2;
					final int alternatives=(int)(rest%16);
					rest/=16;
					final int warnings=(int)(rest%256);
					rest/=256;
					final int claims=(int)(rest%32);
					rest/=32;
					final int source=(int)(rest%SOURCE_LABELS.length);
					rest/=SOURCE_LABELS.length;
					final int margin=(int)(rest%MARGIN_LABELS.length);
					rest/=MARGIN_LABELS.length;
					final int truncation=(int)(rest%4);
					rest/=4;
					final int coverage=(int)(rest%COVERAGE_LABELS.length), query=(int)(rest/COVERAGE_LABELS.length);
					int claimFlags=0;
					if((claims&1)!=0){claimFlags|=SAME_ACCESSION_SOURCE;}
					if((claims&2)!=0){claimFlags|=SAME_FAMILY_SOURCE;}
					if((claims&4)!=0){claimFlags|=OTHER_FAMILY_SOURCE;}
					if((claims&8)!=0){claimFlags|=UNKNOWN_RFAM;}
					if((claims&16)!=0){claimFlags|=RELATED_SOURCE_CLASS;}
					b.clear().append(query==roster.queried() ? "ALL_QUERIED" : roster.query(query).accession).tab();
					b.append(COVERAGE_LABELS[coverage]).tab().append(truncation).tab().append(MARGIN_LABELS[margin]).tab().append(SOURCE_LABELS[source]).tab();
					names(b, claimFlags).tab();
					names(b, warnings<<9).tab();
					for(int bit=1; bit<=8; bit<<=1){b.append((alternatives&bit)!=0).tab();}
					b.append(edge!=0).tab().append(events.get(key)).tab().append(groups.contains(key) ? groups.get(key) : 0);
					out.writer.print(b.append("\tDESCRIPTIVE_BINS_NONADDITIVE_GROUPS\n"));
				}
			}
		}
		final LongIntMap events=new LongIntMap(16), groups=new LongIntMap(16);
		final LongList groupKeys=new LongList();
		int currentGroup=-1;
	}
	private static int coverageBin(final int span, final int length){
		assert(span>0 && span<=length) : "The geometry loader validates model coordinates before display binning";
		if(span==length){return 4;}
		if(span*100L<length*50L){return 0;}
		if(span*100L<length*80L){return 1;}
		if(span*100L<length*95L){return 2;}
		return 3;
	}
	private static int marginBin(final double margin){
		for(int i=0; i<MARGIN_LIMITS.length; i++){if(margin<MARGIN_LIMITS[i]){return i;}}
		return MARGIN_LIMITS.length;
	}

	private static void finishGroup(final long[][] counts, final int[] flags, final boolean[] seen){
		assert(counts.length==flags.length && flags.length==seen.length) : "Per-accession group scratch uses one slot per summary row";
		for(int i=0; i<counts.length; i++){if(seen[i]){count(counts[i], flags[i]);}}
		Arrays.fill(flags, 0);
		Arrays.fill(seen, false);
	}
	private static void count(final long[] counts, final int flags){
		assert(counts.length==NAMES.length+1) : "One total and one counter for every named evidence bit are required";
		counts[0]++;
		for(int i=0; i<NAMES.length; i++){if((flags&(1<<i))!=0){counts[i+1]++;}}
	}
	static ByteBuilder names(final ByteBuilder b, final int flags){
		if(flags==0){return b.append("NONE");}
		boolean first=true;
		for(int i=0; i<NAMES.length; i++){
			if((flags&(1<<i))==0){continue;}
			if(!first){b.append(';');}
			b.append(NAMES[i]);
			first=false;
		}
		return b;
	}
	/** Reads this writer's complete named bit set without allocating field Strings.
	 * Unknown, empty or repeated names fail so a new source warning cannot silently
	 * become an absent warning in a downstream decision. The caller owns the parser.
	 */
	static int readFlags(final LineParser1 parser, final int column){
		if(parser.termEquals("NONE", column)){return 0;}
		parser.setBounds(column);
		final byte[] data=parser.line();
		final int limit=parser.b();
		int result=0;
		for(int from=parser.a(); from<limit; ){
			int to=from;
			while(to<limit && data[to]!=';'){to++;}
			int bit=0;
			for(int i=0; i<NAMES.length && bit==0; i++){
				final String name=NAMES[i];
				if(to-from!=name.length()){continue;}
				boolean equal=true;
				for(int j=0; equal && j<name.length(); j++){equal=data[from+j]==name.charAt(j);}
				if(equal){bit=1<<i;}
			}
			require(bit!=0 && (result&bit)==0, "Unknown, empty or repeated saved source flag");
			result|=bit;
			require(to==limit || to+1<limit, "Trailing empty saved source flag");
			from=to+1;
		}
		require(result!=0, "Empty saved source flag field");
		return result;
	}

	/*--------------------------------------------------------------*/
	/*----------------          Claim parsing       ----------------*/
	/*--------------------------------------------------------------*/

	/** Decode the TSV layer before splitting GFF attributes, then decode each value.
	 * An escaped semicolon in a value never becomes another attribute. Repeated
	 * contradictory booleans retain both observations; free-text Note is not evidence.
	 */
	private static final class Claims{
		int read(final byte[] bytes, final int from, final int to){
			decode(bytes, from, to, raw.clear());
			int flags=0;
			for(int start=0; start<raw.length(); ){
				int end=start;
				while(end<raw.length() && raw.array[end]!=';'){end++;}
				int equal=start;
				while(equal<end && raw.array[equal]!='='){equal++;}
				if(equal<end){
					decode(raw.array, start, equal, key.clear());
					final boolean pseudo=sameText(key, "pseudo"), partial=sameText(key, "partial");
					if(pseudo || partial){
						decode(raw.array, equal+1, end, value.clear());
						if(sameText(value, "true") || sameText(value, "1") || sameText(value, "t")){flags|=pseudo ? PSEUDO_TRUE : PARTIAL_TRUE;}
						else if(sameText(value, "false") || sameText(value, "0") || sameText(value, "f")){flags|=pseudo ? PSEUDO_FALSE : PARTIAL_FALSE;}
						else{flags|=UNKNOWN_BOOLEAN;}
					}else if(sameText(key, "pseudogene")){flags|=PSEUDOGENE_TAG;}
					else if(sameText(key, "exception")){flags|=EXCEPTION;}
					else if(sameText(key, "gene_biotype") || sameText(key, "transcript_biotype") || sameText(key, "biotype")){
						decode(raw.array, equal+1, end, value.clear());
						if(contains(value, "pseudogene")){flags|=PSEUDO_BIOTYPE;}
					}
				}
				start=end+1;
			}
			return flags;
		}
		final ByteBuilder raw=new ByteBuilder(), key=new ByteBuilder(), value=new ByteBuilder();
	}
	private static void decode(final byte[] data, final int from, final int to, final ByteBuilder out){
		assert(from>=0 && from<=to && to<=data.length) : "Percent decoding is confined to one framed field or GFF token";
		for(int i=from; i<to; i++){
			if(data[i]!='%'){out.append(data[i]);}
			else{
				require(i+2<to, "Incomplete percent escape in source claim");
				final int a=hex(data[++i]), b=hex(data[++i]);
				require(a>=0 && b>=0, "Invalid percent escape in source claim");
				out.append((byte)((a<<4)|b));
			}
		}
	}
	private static int hex(final byte value){return value>='0' && value<='9' ? value-'0' : value>='A' && value<='F' ? value-'A'+10 : value>='a' && value<='f' ? value-'a'+10 : -1;}
	private static boolean sameText(final ByteBuilder b, final String text){return b.length()==text.length() && at(b, text, 0);}
	private static boolean contains(final ByteBuilder b, final String text){
		for(int i=0; i+text.length()<=b.length(); i++){if(at(b, text, i)){return true;}}
		return false;
	}
	private static boolean at(final ByteBuilder b, final String text, final int start){
		for(int i=0; i<text.length(); i++){
			final int c=b.array[start+i];
			if((c>='A' && c<='Z' ? c+32 : c)!=text.charAt(i)){return false;}
		}
		return true;
	}
	private static boolean same(final LineParser1 p, final int field, final byte[] bytes){
		p.setBounds(field);
		if(bytes==null || p.b()-p.a()!=bytes.length){return false;}
		for(int i=0; i<bytes.length; i++){if(p.line()[p.a()+i]!=bytes[i]){return false;}}
		return true;
	}
	private static String text(final LineParser1 p, final int field){
		p.setBounds(field);
		return new String(p.line(), p.a(), p.b()-p.a(), StandardCharsets.UTF_8);
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------         Retained state       ----------------*/
	/*--------------------------------------------------------------*/

	static final class SourceInput{
		SourceInput(final InputState state_){state=state_;}
		final InputState state;
		final ObjectMap<String,SourceFile> files=new ObjectMap<String,SourceFile>(String.class, SourceFile.class);
	}
	static final class SourceFile{
		SourceFile(final SourceInput input_, final byte[] path_){
			input=input_;
			path=path_;
		}
		final SourceInput input;
		final byte[] path;
		final LongIntMap lines=new LongIntMap(16);
	}
	static final class SourceRow{
		SourceRow(final SourceFile source_, final long line_, final Target target_, final Target type_, final int index_, final LineParser1 p){
			source=source_;
			line=line_;
			target=target_;
			type=type_;
			index=index_;
			start=p.parseLong(5);
			stop=p.parseLong(6);
			strand=p.termEquals("+", 7) ? '+' : p.termEquals("-", 7) ? '-' : '?';
			assessable=p.termEquals("ASSESSABLE", 12);
			require(!assessable || (start>0 && stop>=start && strand!='?' && p.termEquals("VALID_INTERVAL", 10)),
				"An assessable historical row needs valid coordinates, range and strand");
			selection=p.parseInt(9);
			p.setBounds(4);
			final boolean rna=p.b()-p.a()>=3 && p.line()[p.b()-3]=='R' && p.line()[p.b()-2]=='N' && p.line()[p.b()-1]=='A';
			kind=p.termEquals("STRUCTURAL_CONTEXT", 12) || p.termEquals("gene", 4) || p.termEquals("transcript", 4) || p.termEquals("mRNA", 4) ? CONTEXT
				: rna ? RNA : p.termEquals("exon", 4) || p.termEquals("intron", 4) ? COMPONENT : OTHER;
			parents=p.parseInt(13);
			final int rfam=p.parseInt(14);
			sameExpected=p.parseInt(15);
			oppositeExpected=p.parseInt(16);
			invalidExpected=p.parseInt(17);
			require(parents>=0 && rfam>=0 && sameExpected>=0 && oppositeExpected>=0 && invalidExpected>=0,
				"Historical claim/link counts cannot be negative");
			models=rfam==0 ? NO_MODELS : new Model[rfam];
			if(p.termEquals("STRUCTURAL_CONTEXT", 12)){flags|=STRUCTURAL_CONTEXT;}
			else if(!assessable){flags|=INVALID;}
			if(p.termEquals("true", 11)){flags|=INCOMPATIBLE_ID;}
		}
		final SourceFile source;
		final Target target, type;
		final long line, start, stop;
		final int index, kind, selection, parents, sameExpected, oppositeExpected, invalidExpected;
		final char strand;
		final boolean assessable;
		final Model[] models;
		int flags, modelsSeen, parentsSeen, sameSeen, oppositeSeen, invalidSeen;
		boolean attributesSeen;
	}
	private final NcrnaReferenceGeometry geometry;
	private final NcrnaReferenceProfiles profiles;
	private final Roster roster;
	private final int catalogInputs;
	private final ArrayList<SourceRow> rows=new ArrayList<SourceRow>();
	private final ObjectMap<String,SourceInput> inputs=new ObjectMap<String,SourceInput>(String.class, SourceInput.class);
	private final int[] flags, sameRows, oppositeRows;
	private final int[][] originCounts;
	private final boolean[] rfamSeen;
	private final IntList linkSources=new IntList(), linkHits=new IntList();
	private static final Model[] NO_MODELS=new Model[0];
	private static final int RNA=0, COMPONENT=1, CONTEXT=2, OTHER=3;
	private static final String[] KINDS={"RNA_ROW", "COMPONENT_ROW", "GENE_OR_TRANSCRIPT_CONTEXT", "OTHER_ROW"};
	static final int FULL=1, SHORT=2, TRUNCATED=4, EDGE=8, INVALID=16,
		SAME_ACCESSION_SOURCE=32, SAME_FAMILY_SOURCE=64, OTHER_FAMILY_SOURCE=128, UNKNOWN_RFAM=256,
		PSEUDO_TRUE=512, PSEUDO_FALSE=1024, PSEUDOGENE_TAG=2048, PSEUDO_BIOTYPE=4096,
		PARTIAL_TRUE=8192, PARTIAL_FALSE=16384, EXCEPTION=32768, UNKNOWN_BOOLEAN=65536,
		OPPOSITE_SOURCE=131072, OTHER_MODEL=262144, OTHER_FAMILY_MODEL=524288,
		OPPOSITE_MODEL=1048576, NO_SAME_SOURCE=2097152, STRUCTURAL_CONTEXT=4194304,
		INCOMPATIBLE_ID=8388608, DANGLING_PARENT=16777216, UNTRUNCATED=33554432,
		RELATED_SOURCE_CLASS=67108864, NO_SOURCE_RFAM=134217728, SAME_FAMILY_MODEL=268435456, SAME_MODEL=536870912;
	private static final String[] NAMES={"FULL_MODEL_SPAN", "SHORTER_MODEL_SPAN", "NATIVE_TRUNCATED", "PHYSICAL_CONTIG_EDGE", "UNASSESSABLE_RANGE_OR_SOURCE",
		"SOURCE_SAME_ACCESSION", "SOURCE_SAME_FAMILY_OTHER_ACCESSION", "SOURCE_OTHER_KNOWN_FAMILY", "SOURCE_UNKNOWN_RFAM",
		"SOURCE_PSEUDO_TRUE", "SOURCE_PSEUDO_FALSE", "SOURCE_PSEUDOGENE_ATTRIBUTE", "SOURCE_PSEUDOGENE_BIOTYPE",
		"SOURCE_PARTIAL_TRUE", "SOURCE_PARTIAL_FALSE", "SOURCE_EXCEPTION_ATTRIBUTE", "SOURCE_UNRECOGNIZED_BOOLEAN",
		"OPPOSITE_STRAND_SOURCE_OVERLAP", "OVERLAPPING_OTHER_MODEL", "OVERLAPPING_OTHER_FAMILY_MODEL",
		"OPPOSITE_STRAND_MODEL_OVERLAP", "NO_ASSESSABLE_SAME_STRAND_SOURCE_OVERLAP", "STRUCTURAL_SOURCE_CONTEXT",
		"SOURCE_INCOMPATIBLE_ID", "SOURCE_DANGLING_PARENT", "NATIVE_UNTRUNCATED", "SOURCE_RELATED_CLASS", "NO_OVERLAPPING_SOURCE_RFAM",
		"OVERLAPPING_SAME_FAMILY_OTHER_MODEL", "OVERLAPPING_SAME_MODEL"};
	private static final String[] COVERAGE_LABELS={"0_to_below50_percent", "50_to_below80_percent", "80_to_below95_percent", "95_to_below100_percent", "100_percent"};
	private static final String[] MARGIN_LABELS={"below0", "0_to_below5", "5_to_below10", "10_to_below20", "20_to_below50", "50_to_below100", "100_or_more"};
	private static final double[] MARGIN_LIMITS={0, 5, 10, 20, 50, 100};
	private static final String[] SOURCE_LABELS={"UNASSESSABLE_EVENT", "NO_AVAILABLE_GFF", "AVAILABLE_GFF_NO_ASSESSABLE_SAME_STRAND_OVERLAP", "RNA_ROW_OVERLAP", "NONRNA_ROW_OVERLAP"};
	private static final String ATTRIBUTE_SUMMARY_HEADER="inputs_observed\tcatalog_inputs\tretained_source_rows\tmatched_attribute_rows\tsource_paths_read\tunique_selected_lines\tsource_lines_read\tevidence_status\tstatus";
	static final String[] OUTPUT_SUFFIXES={".event_strata.tsv.gz", ".source_strata.tsv.gz", ".source_links.tsv.gz", ".strata_counts.tsv", ".joint_strata.tsv", ".strata_summary.tsv"};
}
