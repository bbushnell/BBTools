package prok;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import map.ObjectIntMap;
import ml.CellNet;
import ml.CellNetParser;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.NcrnaTextIO.SequenceInput;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.StreamerFactory;
import structures.ByteBuilder;
import structures.ListNum;

/** Conditional endpoint assay on saved matched REAL calls, never synthetic draws.
 * Source ordinals, winning models and physical-end provenance come from saved
 * OFF GFFs. Truth only grades coordinates. One sequence pass supplies features
 * for every seed; each available end is inferred once, then both margins reuse
 * its saved FP32 outputs for single-end comparisons and complete paired seeds.
 * Missing centered models omit that seed/end from both arms, never impute it.
 * Detection, trimming eligibility and path selection are
 * NOT reproduced. CLI/stream ownership follows native sample/learning tools.
 * @author Brian, Ganyu
 */
public final class EvaluateNcrnaRealEndpoints{
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, EvaluateNcrnaRealEndpoints.class, false);
		try{
			final EvaluateNcrnaRealEndpoints tool=new EvaluateNcrnaRealEndpoints(pp.args);
			tool.process();
			pp.outstream.println(tool.geometryOnly ? "NCRNA_REAL_CALL_GEOMETRY_COMPLETE_NO_ALIGNMENT_NO_INFERENCE"
				: tool.gateVectors ? "NCRNA_CENTERED_REAL_GATE_VECTORS_COMPLETE_NO_GATE_FIT"
				: "NCRNA_REAL_ENDPOINT_ASSAY_COMPLETE_CONDITIONAL_NOT_CALLER_REPRODUCTION");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Each invocation binds one entire panel and its declared completed seed/ends. */
	public EvaluateNcrnaRealEndpoints(final String[] args){
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT);
			final String b=pair.length>1 ? pair[1] : null;
			if(a.equals("mode")){
				require("geometry".equals(b) || "evaluate".equals(b) || "gatevectors".equals(b), "mode=geometry, evaluate or gatevectors is required");
				geometryOnly="geometry".equals(b);
				gateVectors="gatevectors".equals(b);
			}else if(a.equals("family")){familyName=b;}
			else if(a.equals("folds")){foldMap=b;}
			else if(a.equals("fold")){heldFold=Integer.parseInt(b);}
			else if(a.equals("assemblyaudit")){assemblyAudit=b;}
			else if(a.equals("sources")){sources=b;}
			else if(a.equals("sourceformat")){sourceFormat=b;}
			else if(a.equals("fastalist")){fastaList=b;}
			else if(a.equals("gff")){gff=b;}
			else if(a.equals("familyroot")){familyRoot=b;}
			else if(a.equals("ref")){ref=b;}
			else if(a.equals("hbm")){hbm=b;}
			else if(a.equals("seeds")){seeds=b;}
			else if(a.equals("modeloffsets")){offsets=b;}
			else if(a.equals("oldtables")){tablePaths[0]=b.split(",");}
			else if(a.equals("newtables")){tablePaths[1]=b.split(",");}
			else if(a.equals("oldnets")){netPaths[0]=b.split(",");}
			else if(a.equals("newnets")){netPaths[1]=b.split(",");}
			else if(a.equals("seedlabels")){labels=b.split(",");}
			else if(a.equals("panel")){panel=b;}
			else if(a.equals("expected")){expected=Parse.parseIntKMG(b);}
			else if(a.equals("witnesses")){expectedWitnesses=Parse.parseIntKMG(b);}
			else if(a.equals("inputs")){expectedInputs=Parse.parseIntKMG(b);}
			else{require(parser.parse(arg, a, b), "Unknown real-endpoint argument: "+arg);}
		}
		in=parser.in1;
		out=parser.out1;
		require(in!=null && out!=null && sources!=null && fastaList!=null && (gateVectors || gff!=null) && familyRoot!=null
			&& ref!=null && hbm!=null && seeds!=null && panel!=null
			&& expected>0 && (gateVectors || expectedWitnesses>=expected) && expectedInputs>0
			&& ("native".equals(sourceFormat) || "legacy".equals(sourceFormat)),
			"Bind witnesses, source roster/format, ordered FASTAs, saved GFF, library and exact counts");
		require("cobalamin_rf00174".equals(familyName) || "cobalamin_rf01482".equals(familyName),
			"This assay supports the explicitly bound main or second cobalamin member");
		member=familyName.equals("cobalamin_rf00174") ? "RF00174" : "RF01482";
		require(!parser.overwrite && !parser.append && parser.in2==null && parser.out2==null && parser.maxReads<0 && parser.samplerate==1,
			"Read the complete declared panel into fresh outputs; sampling and replacement are forbidden");
		if(gateVectors){
			validateGateInputs();
			return;
		}
		require(foldMap==null && heldFold==-2, "folds/fold belong only to the real-call gate adapter");
		if(geometryOnly){
			require(!"legacy".equals(sourceFormat) || assemblyAudit!=null,
				"Legacy geometry requires assemblyaudit= to map exact FASTA paths to canonical fold identities");
			require(offsets==null && labels==null && tablePaths[0]==null && tablePaths[1]==null
				&& netPaths[0]==null && netPaths[1]==null,
				"Geometry mode never consumes learned offsets, tables or networks; do not silently pass unused resources");
			require(Tools.testInputFiles(false, true, in, sources, fastaList, gff, ref, hbm, seeds), "Unreadable real-call geometry binding");
			require(Tools.testOutputFiles(false, false, false, out+".geometry.tsv"), "Geometry output must be fresh");
			return;
		}
		require(assemblyAudit==null, "assemblyaudit= belongs only to geometry mode");
		require(offsets!=null && labels!=null && labels.length>0,
			"Evaluation additionally needs offsets, paired tables/networks and seed labels");
		final HashSet<String> distinct=new HashSet<String>();
		for(final String label:labels){require(!label.isEmpty() && distinct.add(label), "Seed labels must be nonempty and unique");}
		for(int arm=0; arm<2; arm++){
			require(tablePaths[arm]!=null && tablePaths[arm].length==6 && netPaths[arm]!=null && netPaths[arm].length==2*labels.length,
				"Each arm needs six tables in end/probe order and one start/stop net pair per seed label");
			require(Tools.testInputFiles(false, true, tablePaths[arm]), "Unreadable learned table payload");
			for(final String path:netPaths[arm]){
				require(arm==1 && path.equals("none") || Tools.testInputFiles(false, true, path), "Unreadable learned network payload");
			}
		}
		for(int seed=0; seed<labels.length; seed++){
			require(!netPaths[1][seed*2].equals("none") || !netPaths[1][seed*2+1].equals("none"), "Every seed needs at least one completed centered model");
		}
		require(Tools.testInputFiles(false, true, in, sources, fastaList, gff, ref, hbm, seeds, offsets), "Unreadable real-call binding");
		for(final String suffix:SUFFIXES){require(Tools.testOutputFiles(false, false, false, out+suffix), "Output must be fresh: "+out+suffix);}
	}

	private void validateGateInputs(){
		require(gff==null && expectedWitnesses==-1 && foldMap!=null && heldFold>=-1 && heldFold<=4,
			"Gate vectors consume a closed geometry cache, not another witness/GFF join; bind folds= and fold=-1(outer) or0..4");
		require(!"legacy".equals(sourceFormat) || assemblyAudit!=null,
			"Legacy real gate rows need the exact canonical-assembly audit used to create their geometry cache");
		require(offsets!=null && labels!=null && labels.length>0 && tablePaths[0]==null && netPaths[0]==null
			&& tablePaths[1]!=null && tablePaths[1].length==6 && netPaths[1]!=null && netPaths[1].length==2*labels.length,
			"Centered gate extraction uses only six centered tables and one complete first-stage pair per declared seed");
		final HashSet<String> seen=new HashSet<String>();
		for(final String label:labels){require(label.matches("[A-Za-z0-9_-]+") && seen.add(label), "Unique seed labels must also be safe output filename components");}
		require(Tools.testInputFiles(false, true, in, sources, fastaList, ref, hbm, seeds, offsets, foldMap)
			&& Tools.testInputFiles(false, true, tablePaths[1]) && Tools.testInputFiles(false, true, netPaths[1]),
			"Gate inputs require complete readable centered payloads; missing network placeholders are not allowed");
		for(final String label:labels){
			for(final String end:ENDS){
				for(final String suffix:GATE_SUFFIXES){
					require(Tools.testOutputFiles(false, false, false, out+".s"+label+"."+end+suffix), "Every real gate output must be fresh");
				}
			}
		}
		require(Tools.testOutputFiles(false, false, false, out+".features.tsv", out+".selection.tsv"), "Gate features and fold accounting must be fresh");
	}

	/*--------------------------------------------------------------*/
	/*----------------          Processing          ----------------*/
	/*--------------------------------------------------------------*/

	/** Input joins and feature extraction finish before any output is opened. */
	public void process(){
		loadFamily();
		loadSources();
		if(gateVectors){loadGateGeometry();}
		else{
			loadWitnesses();
			bindGff();
		}
		if(geometryOnly){
			for(final Source source:orderedSources){if(!source.calls.isEmpty()){extract(source, null);}}
			for(final Call call:calls){require(call.length>0, "Missing actual FASTA contig for "+call.input+"/"+call.seqid);}
			writeGeometry();
			return;
		}
		loadOffsets();
		final NcrnaPositionalBoundaryFeatures[] features={gateVectors ? null : loadFeatures(0), loadFeatures(1)};
		for(int end=0; !gateVectors && end<2; end++){
			for(int probe=0; probe<3; probe++){
				final RrnaPositionalKmerTable.Table a=geometryTables[0][end][probe];
				final RrnaPositionalKmerTable.Table b=geometryTables[1][end][probe];
				require(a.k==b.k && a.offset==b.offset && a.radius==b.radius && a.alpha==b.alpha && a.records==b.records,
					"Offset comparison must keep table geometry, smoothing and teaching population fixed");
			}
		}
		final boolean upper=Read.TO_UPPER_CASE;
		Read.TO_UPPER_CASE=true;
		try{
			for(final Source source:orderedSources){if(!source.calls.isEmpty()){extract(source, features);}}
		}finally{Read.TO_UPPER_CASE=upper;}
		for(final Call call:calls){require(call.features!=null, "Missing actual FASTA contig for "+call.input+"/"+call.seqid);}
		writeFeatures();
		if(gateVectors){writeGateVectors();}else{evaluate();}
	}

	private void loadFamily(){
		final NcrnaFamilyLoader.Options options=new NcrnaFamilyLoader.Options();
		options.resources.root=familyRoot;
		options.resources.consensus=ref;
		options.resources.models=hbm;
		options.resources.seeds=seeds;
		options.boundaryOverride=false;
		options.endpoints=false;
		final FamilyResourceResolver resolver=new FamilyResourceResolver();
		resolver.bindRoot(familyName, familyRoot);
		family=NcrnaFamilyLoader.load(familyName, options, resolver);
		require(family.minLen==109 && family.maxLen==Integer.MAX_VALUE && family.physicalEdgeMinLen==0
			&& family.prok16sPartialMinLen==0 && (family.modelThresholds==null || !family.modelThresholds.hasLengths()),
			"This conditional cobalamin assay binds the ordinary 109-base minimum and no clipped/model length override");
		modelIndex=new ObjectIntMap<String>(family.modelNames.length*2, String.class);
		plateaus=new float[family.library.length];
		for(int i=0; i<family.modelNames.length; i++){
			require(!modelIndex.contains(family.modelNames[i]) && family.models[i]!=null
				&& Arrays.equals(family.library[i], family.models[i].original), "Winning model IDs must be unique and HBMs exactly paired");
			modelIndex.put(family.modelNames[i], i);
			plateaus[i]=TrnaBoundaryFeatures.modelMaxCoverage(family.models[i]);
			require(plateaus[i]>0, "The actual selected HBM needs positive training coverage");
		}
	}

	private void loadSources(){
		final HashMap<String, Source> byFasta=new HashMap<String, Source>();
		final LineParser1 p=new LineParser1('\t');
		final boolean nativeFormat="native".equals(sourceFormat);
		try(Input input=new Input(sources)){
			RrnaResourceIO.header(input.reader, nativeFormat ? SOURCE_HEADER : LEGACY_HEADER, sources);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				require(p.terms()==(nativeFormat ? 5 : 11), "Source roster row width changed");
				final Source source=new Source(p.parseString(nativeFormat ? 2 : 1), p.parseString(nativeFormat ? 3 : 7),
					p.parseString(nativeFormat ? 1 : 2));
				source.savedOrdinal=p.parseInt(0);
				require(byFasta.put(source.fasta, source)==null && byInput.put(source.id, source)==null, "Source paths and input IDs must be one-to-one");
			}
		}
		if(assemblyAudit!=null){
			try(Input input=new Input(assemblyAudit)){
				RrnaResourceIO.header(input.reader, NcrnaEndpointFoldInputs.AUDIT_LABEL, assemblyAudit);
				RrnaResourceIO.header(input.reader, NcrnaEndpointFoldInputs.AUDIT_HEADER, assemblyAudit);
				for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
					p.set(row);
					require(p.terms()==5 && p.termEquals("EXCLUDED_ASSEMBLY", 4) && p.length(3)>0
						&& !p.termEquals(".", 3) && !p.termEquals("NA", 3), "Canonical assembly audit needs exact accepted bindings");
					final Source source=byFasta.get(p.parseString(2));
					if(source==null){continue;}
					final String assembly=p.parseString(3);
					require(!source.canonicalAssembly || source.assembly.equals(assembly), "One FASTA cannot bind to multiple canonical assemblies");
					require(!nativeFormat || source.assembly.equals(assembly), "Native canonical identity must agree with the audit");
					source.assembly=assembly;
					source.canonicalAssembly=true;
				}
			}
		}
		try(Input input=new Input(fastaList)){
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				require(p.terms()==1 && p.length(0)>0, "Ordered caller FASTA list contains one literal path per row");
				final Source source=byFasta.get(p.parseString(0));
				require(source!=null && source.ordinal<0, "Every caller FASTA needs one unique explicit roster binding");
				require(assemblyAudit==null || source.canonicalAssembly, "Every selected legacy FASTA must resolve through the canonical assembly audit");
				source.ordinal=orderedSources.size();
				require(!nativeFormat || source.savedOrdinal==source.ordinal, "Native source ordinal must equal saved caller input order");
				orderedSources.add(source);
			}
		}
		require(orderedSources.size()==expectedInputs, "Caller input list did not close at the accepted panel size");
	}

	private void loadWitnesses(){
		final HashSet<String> references=new HashSet<String>();
		final LineParser1 p=new LineParser1('\t');
		int rows=0;
		try(Input input=new Input(in)){
			RrnaResourceIO.header(input.reader, WITNESS_HEADER, in);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				require(p.terms()==12, "Saved matched witness schema changed");
				rows++;
				if(!p.termEquals(member, 2)){continue;}
				final Call call=new Call();
				call.input=p.parseString(0);
				call.reference=p.parseString(1);
				call.modelName=p.parseString(3);
				call.seqid=p.parseString(4);
				call.low=p.parseInt(5);
				call.high=p.parseInt(6);
				call.truthLow=p.parseInt(7);
				call.truthHigh=p.parseInt(8);
				require(p.termEquals("+", 9) || p.termEquals("-", 9), "Real calls need an explicit genomic strand");
				call.minus=p.termEquals("-", 9);
				final Source source=byInput.get(call.input);
				require(source!=null && source.ordinal>=0 && modelIndex.contains(call.modelName), "Witness source and exact winning model must exist in the frozen caller bindings");
				call.model=modelIndex.get(call.modelName);
				call.ordinal=source.ordinal;
				call.assembly=source.assembly;
				call.row=calls.size()+1;
				require(call.low>=1 && call.high>=call.low && call.truthLow>=1 && call.truthHigh>=call.truthLow,
					"Witness coordinates are one-based inclusive genomic spans");
				require(p.parseLong(10)==(call.minus ? (long)call.truthHigh-call.high : (long)call.low-call.truthLow)
					&& p.parseLong(11)==(call.minus ? (long)call.truthLow-call.low : (long)call.high-call.truthHigh),
					"Witness signed errors must equal caller-minus-truth in transcript orientation");
				require(references.add(call.input+'\t'+call.reference) && source.keys.put(call.key(), call)==null,
					"Saved matching must consume each source-qualified truth and physical call once");
				ArrayList<Call> contig=source.calls.get(call.seqid);
				if(contig==null){
					contig=new ArrayList<Call>();
					source.calls.put(call.seqid, contig);
				}
				contig.add(call);
				calls.add(call);
			}
		}
		require(rows==expectedWitnesses && calls.size()==expected, "Full witness and selected-member denominators must match the accepted original cohorts");
	}

	/** Read the already-qualified source/geometry cache. Fold assignment is
	 * assembly-level and explicit even for assemblies with no teaching locus.
	 * All cache records are validated before selecting; no unknown assembly can
	 * silently default to training. Truth supplies labels only after extraction. */
	private void loadGateGeometry(){
		assert(gateVectors) : "Only the centered gate adapter consumes cached geometry";
		final ObjectIntMap<String> assigned=new ObjectIntMap<String>(String.class);
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(foldMap)){
			RrnaResourceIO.header(input.reader, NcrnaEndpointFoldInputs.ASSEMBLY_HEADER, foldMap);
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==4 && p.length(0)>0 && p.termStartsWith("/", 1), "Assembly folds need exact identities and an absolute canonical FASTA");
				final String assembly=p.parseString(0);
				final int fold=p.parseInt(3);
				require(!assigned.contains(assembly) && (p.termEquals("OUTER_TEST", 2) ? fold==-1
					: p.termEquals("INTERNAL_FOLD", 2) && fold>=0 && fold<5), "One assembly has one explicit outer/internal fold role");
				assigned.put(assembly, fold+1);
			}
		}
		require(assigned.size()>0, "Real gate extraction requires a nonempty complete assembly mapping");
		final HashSet<String> references=new HashSet<String>();
		int rows=0;
		boolean complete=false;
		try(Input input=new Input(in)){
			final boolean legacyMaximum=readGeometryHeader(input.reader, false, in);
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				if(p.termEquals("#COMPLETE", 0)){
					require(p.terms()==2 && p.parseInt(1)==expected && rows==expected && input.reader.nextLine()==null,
						"Geometry cache must end once at the accepted full member count");
					complete=true;
					break;
				}
				require(p.terms()==16 && p.parseInt(0)==++rows && p.termEquals(member, 4), "Gate cache rows preserve dense original identity and explicit member");
				final Call call=new Call();
				call.row=rows;
				call.input=p.parseString(1);
				call.assembly=p.parseString(2);
				call.reference=p.parseString(3);
				call.modelName=p.parseString(5);
				call.seqid=p.parseString(6);
				require(p.termEquals("+", 7) || p.termEquals("-", 7), "Cached real calls retain an explicit strand");
				call.minus=p.termEquals("-", 7);
				call.expectedLength=p.parseInt(8);
				call.fixed=p.parseInt(9);
				final int start=p.parseInt(10), stop=p.parseInt(11);
				final int truthStart=p.parseInt(12), truthStop=p.parseInt(13);
				require(call.expectedLength>0 && call.fixed>=0 && call.fixed<=3 && start>=0 && stop>=start && stop<call.expectedLength
					&& truthStart>=0 && truthStop>=truthStart && truthStop<call.expectedLength
					&& p.parseInt(14)==family.minLen && NcrnaTextIO.maximumLength(p, 15, legacyMaximum)==family.maxLen,
					"Cached pairs must fit actual contig bounds and retain ordinary family length guards; raw short spans remain in the cohort");
				require(((call.fixed&1)==0 || start==0) && ((call.fixed&2)==0 || stop==call.expectedLength-1),
					"Cached fixed bits preserve physical ends, not inferred edge-contact eligibility");
				final Source source=byInput.get(call.input);
				require(source!=null && source.ordinal>=0 && source.assembly.equals(call.assembly) && modelIndex.contains(call.modelName)
					&& assigned.contains(call.assembly), "Every cache row needs its exact original FASTA/source, canonical assembly and winning model");
				call.model=modelIndex.get(call.modelName);
				call.ordinal=source.ordinal;
				call.low=call.minus ? call.expectedLength-stop : start+1;
				call.high=call.minus ? call.expectedLength-start : stop+1;
				call.truthLow=call.minus ? call.expectedLength-truthStop : truthStart+1;
				call.truthHigh=call.minus ? call.expectedLength-truthStart : truthStop+1;
				require(references.add(call.input+'\t'+call.reference) && source.keys.put(call.key(), call)==null,
					"Cached source-qualified truth and physical calls must remain one-to-one before fold selection");
				final int fold=assigned.get(call.assembly)-1;
				foldCounts[fold+1]++;
				if(fold!=heldFold){continue;}
				ArrayList<Call> contig=source.calls.get(call.seqid);
				if(contig==null){
					contig=new ArrayList<Call>();
					source.calls.put(call.seqid, contig);
				}
				contig.add(call);
				calls.add(call);
			}
		}
		require(complete && calls.size()==foldCounts[heldFold+1], "Fold selection must conserve every declared held-call record, including an empty chunk");
		require(heldFold!=-1 || calls.size()==expected, "An outer-panel gate invocation cannot silently discard internal calls");
	}

	/** No edge-contact heuristic: only an exact saved GFF match supplies fixed bits. */
	private void bindGff(){
		final LineParser1 p=new LineParser1('\t');
		final boolean[] seen=new boolean[expectedInputs];
		Source source=null;
		try(Input input=new Input(gff)){
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				if(p.termEquals("##input_ordinal0", 0)){
					require(p.terms()==2, "Caller source marker requires its ordinal");
					final int ordinal=p.parseInt(1);
					require(ordinal>=0 && ordinal<seen.length && !seen[ordinal], "Saved caller input markers must be unique and within the declared panel");
					seen[ordinal]=true;
					source=orderedSources.get(ordinal);
					continue;
				}
				if(row.length==0 || row[0]=='#'){continue;}
				require(source!=null && p.terms()==9, "Every GFF row needs an explicit caller input ordinal and nine columns");
				if(!p.termEquals("RNA", 2) || source.keys.isEmpty()){continue;}
				final String key=p.parseString(0)+'\t'+p.parseString(6)+'\t'+p.parseInt(3)+'\t'+p.parseInt(4);
				final Call call=source.keys.get(key);
				if(call==null){continue;}
				final String attributes=p.parseString(8);
				require(call.fixed<0 && call.modelName.equals(attribute(attributes, "model")) && familyName.equals(attribute(attributes, "family")),
					"Saved GFF call must uniquely match witness coordinates, strand, family and full winning model");
				call.fixed=fixedMask(attributes, call.minus);
			}
		}
		for(final boolean marker:seen){require(marker, "Missing caller source marker in saved GFF");}
		for(final Call call:calls){require(call.fixed>=0, "Witness has no exact saved GFF provenance: "+call.input+"/"+call.reference);}
	}

	/** Parses only whole comma/semicolon-delimited GFF note keys; duplicates fail. */
	public static String attribute(final String text, final String key){
		String value=null;
		for(int from=0; from<text.length();){
			int to=from;
			while(to<text.length() && text.charAt(to)!=',' && text.charAt(to)!=';'){to++;}
			if(text.startsWith(key, from) && from+key.length()<to
				&& (text.charAt(from+key.length())==':' || text.charAt(from+key.length())=='=')){
				require(value==null, "Duplicate saved GFF attribute: "+key);
				value=text.substring(from+key.length()+1, to);
			}
			from=to+1;
		}
		return value;
	}

	public static int fixedMask(final String attributes, final boolean minus){
		final String partial=attribute(attributes, "partial");
		if(partial==null){return 0;}
		require(partial.equals("10") || partial.equals("01") || partial.equals("11"), "Orf GFF partial provenance must be 10, 01 or 11");
		final int genomic=(partial.charAt(0)=='1' ? 1 : 0)|(partial.charAt(1)=='1' ? 2 : 0);
		return minus ? ((genomic&1)<<1)|((genomic&2)>>1) : genomic;
	}

	private void loadOffsets(){
		median5=new int[family.modelNames.length];
		median3=new int[median5.length];
		support=new long[median5.length];
		final LineParser1 p=new LineParser1('\t');
		int model=0;
		try(Input input=new Input(offsets)){
			RrnaResourceIO.header(input.reader, NcrnaFixedEndpointTrials.OFFSET_HEADER, offsets);
			for(byte[] row=input.reader.nextLine(); row!=null; row=input.reader.nextLine()){
				p.set(row);
				require(p.terms()==6, "Offset resource needs six explicit fields");
				if(!p.termEquals(member, 0)){continue;}
				require(model<median5.length && p.termEquals(NcrnaSupport.recordName(family.modelNames[model]), 1)
					&& p.parseInt(3)==10 && p.parseLong(2)>=0, "Training offsets must preserve the frozen library order and support floor ten");
				support[model]=p.parseLong(2);
				median5[model]=Math.toIntExact(p.parseLong(4));
				median3[model]=Math.toIntExact(p.parseLong(5));
				require(support[model]>=10 || median5[model]==0 && median3[model]==0, "Unsupported models cannot acquire an offset");
				model++;
			}
		}
		require(model==median5.length, "Offsets must cover every model, including unsupported zero-offset models");
	}

	private NcrnaPositionalBoundaryFeatures loadFeatures(final int arm){
		final NcrnaWeightedEndpointTable[][] weighted=new NcrnaWeightedEndpointTable[2][3];
		final RrnaPositionalKmerTable.Table[][] tables=new RrnaPositionalKmerTable.Table[2][3];
		for(int end=0; end<2; end++){
			for(int probe=0; probe<3; probe++){
				weighted[end][probe]=NcrnaWeightedEndpointTable.loadScores(tablePaths[arm][3*end+probe], member);
				tables[end][probe]=weighted[end][probe].geometry;
			}
		}
		final NcrnaPositionalBoundaryFeatures features=new NcrnaPositionalBoundaryFeatures(tables[0], tables[1], weighted[0], weighted[1], false);
		geometryTables[arm]=tables;
		require(features.radius==10 && features.values[0].length==69, "Retain the qualified radius-ten 69-input feature layout");
		return features;
	}

	private void extract(final Source source, final NcrnaPositionalBoundaryFeatures[] features){
		final FileFormat format=FileFormat.testInput(source.fasta, FileFormat.FASTA, null, true, true);
		try(SequenceInput input=new SequenceInput(StreamerFactory.makeStreamer(format, null, true, -1, false, true, 1))){
			input.stream.start();
			for(ListNum<Read> batch=input.stream.nextList(); batch!=null; batch=input.stream.nextList()){
				for(final Read read:batch){
					final ArrayList<Call> selected=source.calls.get(NcrnaSupport.recordName(read.id));
					if(selected==null){continue;}
					require(read.bases!=null && read.mate==null, "Selected real calls need intact unpaired genomic contigs");
					final float gc=geometryOnly ? 0 : GeneCaller.calcContigGC(read.bases);
					byte[] reverse=null;
					for(final Call call:selected){
						require(call.length==0 && call.high<=read.length() && call.truthHigh<=read.length(), "Repeated contig ID or witness beyond actual FASTA length");
						require(!gateVectors || read.length()==call.expectedLength, "Actual FASTA length changed since the accepted geometry cache");
						if(!geometryOnly && call.minus && reverse==null){
							reverse=read.bases.clone();
							AminoAcid.reverseComplementBasesInPlace(reverse);
						}
						final byte[] bases=!geometryOnly && call.minus ? reverse : read.bases;
						call.length=bases.length;
						call.start=senseStart(call.low, call.high, bases.length, call.minus);
						call.stop=senseStop(call.low, call.high, bases.length, call.minus);
						call.truthStart=senseStart(call.truthLow, call.truthHigh, bases.length, call.minus);
						call.truthStop=senseStop(call.truthLow, call.truthHigh, bases.length, call.minus);
						require(((call.fixed&1)==0 || call.start==0) && ((call.fixed&2)==0 || call.stop==bases.length-1),
							"Saved partial bits must identify physical contig ends in transcript orientation");
						if(geometryOnly){continue;}
						call.proposedStart=(long)call.start-((call.fixed&1)==0 ? median5[call.model] : 0);
						call.proposedStop=(long)call.stop-((call.fixed&2)==0 ? median3[call.model] : 0);
						call.centerBad=geometry(call.proposedStart, call.proposedStop, bases.length, family.minLen, family.maxLen);
						call.centerStart=call.centerBad==0 ? (int)call.proposedStart : call.start;
						call.centerStop=call.centerBad==0 ? (int)call.proposedStop : call.stop;
						call.features=new float[2][2][];
						for(int arm=gateVectors ? 1 : 0; arm<2; arm++){
							if(arm==1 && call.centerBad!=0){continue;}
							call.usable[arm]=features[arm].fill(bases, arm==0 ? call.start : call.centerStart,
								arm==0 ? call.stop : call.centerStop, family.library[call.model], family.models[call.model], plateaus[call.model], gc, call.fixed);
							for(int end=0; end<2; end++){
								if((call.usable[arm]&(1<<end))!=0){call.features[arm][end]=features[arm].values[end].clone();}
							}
						}
						if(gateVectors && call.usable[1]!=0){call.centerGC=NcrnaEndpointGateFeatures.gc(bases, call.centerStart, call.centerStop);}
					}
				}
			}
		}
	}

	public static int senseStart(final int low, final int high, final int length, final boolean minus){
		require(low>=1 && high>=low && high<=length, "Genomic bounds must fit their actual contig");
		return minus ? length-high : low-1;
	}
	public static int senseStop(final int low, final int high, final int length, final boolean minus){
		require(low>=1 && high>=low && high<=length, "Genomic bounds must fit their actual contig");
		return minus ? length-low : high-1;
	}
	/** Mirrors the native joint guard, distinguishing its failure causes. */
	public static int geometry(final long start, final long stop, final int length, final int minimum, final int maximum){
		require(length>0 && minimum>0 && maximum>=minimum, "Joint geometry needs real contig length and caller length bounds");
		return stop<start ? 1 : start<0 || stop>=length ? 2 : stop-start+1<minimum ? 3 : stop-start+1>maximum ? 4 : 0;
	}

	/** Every matched call stays in the cache, including raw spans below the
	 * family minimum. Later policy guards must hold invalid proposals, not remove
	 * difficult rows from the denominator. Truth is never used to derive geometry.
	 */
	private void writeGeometry(){
		assert(geometryOnly && calls.size()==expected) : "Geometry output requires the complete bound matched cohort";
		final ByteBuilder row=new ByteBuilder();
		try(Output output=new Output(out+".geometry.tsv")){
			output.writer.println(GEOMETRY_HEADER);
			for(final Call call:calls){
				row.clear().append(call.row).tab().append(call.input).tab().append(call.assembly)
					.tab().append(call.reference).tab().append(member).tab().append(call.modelName).tab().append(call.seqid)
					.tab().append(call.minus ? "-" : "+").tab().append(call.length).tab().append(call.fixed)
					.tab().append(call.start).tab().append(call.stop).tab().append(call.truthStart).tab().append(call.truthStop)
					.tab().append(family.minLen).tab();
				NcrnaTextIO.appendMaximumLength(row, family.maxLen);
				output.writer.print(row.nl());
			}
			output.writer.print(row.clear().append("#COMPLETE\t").append(calls.size()).nl());
		}
	}

	private void writeFeatures(){
		final ByteBuilder b=new ByteBuilder();
		try(Output output=new Output(out+".features.tsv")){
			output.writer.println("row1\tarm\tend\tstatus\tfeature_bits_0_through_68");
			for(final Call call:calls){
				for(int arm=gateVectors ? 1 : 0; arm<2; arm++){
					for(int end=0; end<2; end++){
						b.clear().append(call.row).tab().append(arm==0 ? "raw" : "center").tab().append(ENDS[end]).tab().append(call.status(arm, end));
						final float[] values=call.features[arm][end];
						if(values!=null){for(final float value:values){b.tab().append(Float.floatToRawIntBits(value));}}
						output.writer.print(b.nl());
					}
				}
			}
			output.writer.print(b.clear().append("#COMPLETE\t").append(calls.size()).nl());
		}
	}

	/** Each real held-call endpoint is scored once. Gate rows are only proposals
	 * that differ from RAW; offset-only moves remain candidates. The audit keeps
	 * deterministic holds and physically invalid proposals in the original cohort.
	 * Later cutoff replay applies single/joint physical guards without inference. */
	private void writeGateVectors(){
		assert(gateVectors && calls.size()==foldCounts[heldFold+1]) : "Gate output requires the exact assembly-held cohort";
		final boolean bgzf=ReadWrite.USE_BGZF;
		final boolean force=ReadWrite.FORCE_BGZIP;
		ReadWrite.USE_BGZF=true;
		ReadWrite.FORCE_BGZIP=true;
		try{
			final float[] referenceGC=new float[family.library.length];
			for(int model=0; model<referenceGC.length; model++){
				referenceGC[model]=NcrnaEndpointGateFeatures.gc(family.library[model], 0, family.library[model].length-1);
			}
			for(int seed=0; seed<labels.length; seed++){
				for(int end=0; end<2; end++){writeGateEnd(seed, end, referenceGC);}
			}
			final ByteBuilder row=new ByteBuilder();
			try(Output output=new Output(out+".selection.tsv")){
				output.writer.println("panel\tmember\trequested_fold\tcache_fold\tcalls\tselected\tschema\tqualification");
				for(int fold=-1; fold<5; fold++){
					output.writer.print(row.clear().append(panel).tab().append(member).tab().append(heldFold).tab().append(fold)
						.tab().append(foldCounts[fold+1]).tab().append(fold==heldFold ? 1 : 0).tab().append(NcrnaCenteredGateFeatures.SCHEMA)
						.append("\tCONDITIONAL_REAL_CALLS_FOLD4_NOT_NESTED_INDEPENDENT\n"));
				}
			}
		}finally{
			ReadWrite.USE_BGZF=bgzf;
			ReadWrite.FORCE_BGZIP=force;
		}
	}

	private void writeGateEnd(final int seed, final int end, final float[] referenceGC){
		final String prefix=out+".s"+labels[seed]+"."+ENDS[end];
		final CellNet net=CellNetParser.load(netPaths[1][2*seed+end], false);
		RrnaEndpointCallerFeatures.Resources.validateNet(net, 69, 22);
		final ByteBuilder row=new ByteBuilder();
		final float[] scores=new float[22], scratch=new float[32], values=new float[35];
		long proposals=0, positive=0, inferred=0;
		try(Output matrix=new Output(prefix+".vectors.tsv.gz"); Output audit=new Output(prefix+".rows.tsv.gz");
				Output saved=new Output(prefix+".scores.tsv.gz"); Output summary=new Output(prefix+".summary.tsv")){
			matrix.writer.println("#dims\t35\t1\t0");
			gateHeader(audit, seed, end);
			audit.writer.println(GATE_HEADER);
			gateHeader(saved, seed, end);
			row.clear().append("row1\tstatus");
			for(int i=0; i<22; i++){row.append("\tscore_bits_").append(i);}
			saved.writer.print(row.nl());
			for(final Call call:calls){
				final float[] features=call.features[1][end];
				final int raw=end==0 ? call.start : call.stop;
				final int center=end==0 ? call.centerStart : call.centerStop;
				long proposal=raw, gateRow=0;
				int neural=0;
				boolean better=false;
				String status=call.status(1, end);
				if(features!=null){
					require(call.centerGC>=0 && call.centerGC<=1 && Float.isFinite(call.centerGC),
						"A vectorized centered span requires its actually measured GC fraction");
					net.applyInput(features);
					net.feedForward();
					for(int i=0; i<22; i++){
						scores[i]=net.getOutput(i);
						require(Float.isFinite(scores[i]), "Cached real first-stage scores cannot contain nonfinite evidence");
					}
					inferred++;
					neural=NcrnaCenteredGateFeatures.adjustment(scores);
					proposal=NcrnaCenteredGateFeatures.fill(scores, features, call.centerGC, referenceGC[call.model],
						family.library[call.model].length, call.start, call.stop, call.centerStart, call.centerStop, end, scratch, values);
					status=proposal==raw ? "SAME_AS_RAW_HOLD" : "PROPOSAL";
					if(proposal!=raw){
						better=NcrnaCenteredGateFeatures.improves(raw, proposal, end==0 ? call.truthStart : call.truthStop);
						gateRow=++proposals;
						if(better){positive++;}
						row.clear();
						for(final float value:values){row.append(value, 18).tab();}
						matrix.writer.print(row.append(better ? 1 : 0).nl());
					}
				}
				row.clear().append(call.row).tab().append(status);
				if(features!=null){for(final float value:scores){row.tab().append(Float.floatToRawIntBits(value));}}
				saved.writer.print(row.nl());
				final int bad=geometry(end==0 ? proposal : call.start, end==1 ? proposal : call.stop,
					call.length, family.minLen, family.maxLen);
				geometryRow(row.clear(), call);
				audit.writer.print(row.tab().append(call.centerStart).tab().append(call.centerStop).tab().append(support[call.model])
					.tab().append(GEOMETRY[call.centerBad]).tab().append(ENDS[end]).tab().append(status).tab().append(gateRow)
					.tab().append(proposal).tab().append(neural).tab().append((long)center-raw).tab().append(better ? 1 : 0)
					.tab().append(GEOMETRY[bad]).nl());
			}
			audit.writer.print(row.clear().append("#COMPLETE\t").append(calls.size()).tab().append(proposals).nl());
			saved.writer.print(row.clear().append("#COMPLETE\t").append(calls.size()).tab().append(inferred).nl());
			summary.writer.println("panel\tmember\tfold\tseed\tend\tcache_calls\theld_calls\tinferred\tproposals\tpositives\tdeterministic_holds\tinputs\tschema\tfeature_order");
			summary.writer.print(row.clear().append(panel).tab().append(member).tab().append(heldFold).tab().append(labels[seed])
				.tab().append(ENDS[end]).tab().append(expected).tab().append(calls.size()).tab().append(inferred).tab().append(proposals)
				.tab().append(positive).tab().append(calls.size()-proposals).tab().append(35).tab().append(NcrnaCenteredGateFeatures.SCHEMA)
				.tab().append(NcrnaCenteredGateFeatures.FEATURE_ORDER).nl());
		}
	}

	private void gateHeader(final Output output, final int seed, final int end){
		assert(gateVectors && heldFold>=-1) : "Real gate provenance must bind an explicit fold or outer role";
		output.writer.print(new ByteBuilder("#NCRNA_CENTERED_REAL_GATE_V1\n#schema\t").append(NcrnaCenteredGateFeatures.SCHEMA)
			.append("\n#panel\t").append(panel).append("\n#member\t").append(member).append("\n#fold\t").append(heldFold)
			.append("\n#geometry\t").append(in).append("\n#folds\t").append(foldMap).append("\n#offsets\t").append(offsets)
			.append("\n#network\t").append(netPaths[1][2*seed+end]).append("\n#seed\t").append(labels[seed])
			.append("\n#end\t").append(ENDS[end]).nl());
	}

	private void geometryRow(final ByteBuilder row, final Call call){
		assert(call.length>0) : "A real-call row must be re-bound to its actual FASTA before output";
		row.append(call.row).tab().append(call.input).tab().append(call.assembly).tab().append(call.reference)
			.tab().append(member).tab().append(call.modelName).tab().append(call.seqid).tab().append(call.minus ? "-" : "+")
			.tab().append(call.length).tab().append(call.fixed).tab().append(call.start).tab().append(call.stop)
			.tab().append(call.truthStart).tab().append(call.truthStop).tab().append(family.minLen).tab();
		NcrnaTextIO.appendMaximumLength(row, family.maxLen);
	}

	/** Recognizes the declared old numeric-maximum schema without weakening any
	 * other header field. Current schemas name the maximum's conditional presence. */
	public static boolean readGeometryHeader(final ByteFile input, final boolean gate, final String path){
		final byte[] header=input.nextLine();
		final boolean legacy=Arrays.equals(header, (gate ? LEGACY_GATE_HEADER : LEGACY_GEOMETRY_HEADER).getBytes(StandardCharsets.US_ASCII));
		require(legacy || Arrays.equals(header, (gate ? GATE_HEADER : GEOMETRY_HEADER).getBytes(StandardCharsets.US_ASCII)),
			"Unexpected real-call geometry header: "+path);
		return legacy;
	}

	private void evaluate(){
		final Stats[][][][][] stats=new Stats[labels.length][2][2][4][2];
		final ByteBuilder b=new ByteBuilder();
		final float[][][] scoreBuffers=new float[2][2][22];
		final float[][][] output=new float[2][2][];
		final int[][] proposed=new int[4][2];
		final int[][] finalEnds=new int[4][2];
		final int[] bad=new int[4];
		try(Output scores=new Output(out+".scores.tsv"); Output trials=new Output(out+".trials.tsv")){
			scores.writer.println("row1\tseed\tarm\tend\tstatus\tscore_bits_0_through_21");
			trials.writer.println("row1\tinput_id\treference_id\tmodel\tseqid\tstrand\tcontig_length\tfixed_mask\ttraining_n\traw_start0\traw_stop0\ttrue_start0\ttrue_stop0\tproposed_center_start0\tproposed_center_stop0\tcenter_geometry\tseed\tmargin\tarm\tproposed_start0\tproposed_stop0\tjoint_geometry\tfinal_start0\tfinal_stop0\td5\td3\tscope\ttarget_end");
			for(int seed=0; seed<labels.length; seed++){
				final CellNet[][] nets=new CellNet[2][2];
				for(int arm=0; arm<2; arm++){
					for(int end=0; end<2; end++){
						if(!available(seed, end)){continue;}
						nets[arm][end]=CellNetParser.load(netPaths[arm][seed*2+end], false);
						RrnaEndpointCallerFeatures.Resources.validateNet(nets[arm][end], 69, 22);
					}
				}
				for(int scope=0; scope<2; scope++){
					for(int margin=0; margin<2; margin++){
						for(int arm=0; arm<4; arm++){
							for(int end=0; end<2; end++){
								if(available(seed, end) && (scope==0 || available(seed, 1-end))){stats[seed][scope][margin][arm][end]=new Stats();}
							}
						}
					}
				}
				for(final Call call:calls){
					for(int arm=0; arm<2; arm++){
						for(int end=0; end<2; end++){
							final float[] values=call.features[arm][end];
							output[arm][end]=null;
							if(values!=null && available(seed, end)){
								final CellNet net=nets[arm][end];
								net.applyInput(values);
								net.feedForward();
								output[arm][end]=scoreBuffers[arm][end];
								for(int i=0; i<22; i++){
									final float value=net.getOutput(i);
									require(Float.isFinite(value), "Real-call network outputs must be finite");
									output[arm][end][i]=value;
								}
							}
							b.clear().append(call.row).tab().append(labels[seed]).tab().append(arm==0 ? "raw" : "center").tab().append(ENDS[end])
								.tab().append(available(seed, end) ? call.status(arm, end) : "MODEL_UNAVAILABLE");
							if(output[arm][end]!=null){for(final float score:output[arm][end]){b.tab().append(Float.floatToRawIntBits(score));}}
							scores.writer.print(b.nl());
						}
					}
					for(int scope=0; scope<2; scope++){
						for(int target=0; target<(scope==0 ? 2 : 1); target++){
							if(!available(seed, target) || scope==1 && !available(seed, 1)){continue;}
							for(int margin=0; margin<2; margin++){
								proposed[0][0]=call.start;
								proposed[0][1]=call.stop;
								proposal(call.start, call.stop, output[0], MARGINS[margin], proposed[1]);
								proposed[2][0]=call.centerStart;
								proposed[2][1]=call.centerStop;
								proposal(call.centerStart, call.centerStop, output[1], MARGINS[margin], proposed[3]);
								if(scope==0){
									final int other=1-target;
									proposed[1][other]=proposed[0][other];
									proposed[3][other]=proposed[2][other];
								}
								for(int arm=0; arm<4; arm++){
									bad[arm]=geometry(proposed[arm][0], proposed[arm][1], call.length, family.minLen, family.maxLen);
									final int baseStart=arm>=2 ? call.centerStart : call.start;
									final int baseStop=arm>=2 ? call.centerStop : call.stop;
									finalEnds[arm][0]=bad[arm]==0 ? proposed[arm][0] : baseStart;
									finalEnds[arm][1]=bad[arm]==0 ? proposed[arm][1] : baseStop;
									writeTrial(trials, b, call, seed, margin, arm, proposed[arm], bad[arm], finalEnds[arm], scope, target);
								}
								for(int arm=0; arm<4; arm++){
									final int comparator=arm==3 ? 1 : 0;
									for(int end=0; end<2; end++){
										if(scope==0 && end!=target){continue;}
										final long truth=end==0 ? call.truthStart : call.truthStop;
										stats[seed][scope][margin][arm][end].add(finalEnds[arm][end]-truth, finalEnds[0][end]-truth,
											finalEnds[comparator][end]-truth, finalEnds[arm][1-end]-(end==0 ? call.truthStop : call.truthStart),
											bad[arm], call, arm, end);
									}
								}
							}
						}
					}
				}
			}
			scores.writer.print(b.clear().append("#COMPLETE\t").append(calls.size()).tab().append(labels.length).nl());
			trials.writer.print(b.clear().append("#COMPLETE\t").append(calls.size()).tab().append(labels.length).nl());
		}
		try(Output summary=new Output(out+".summary.tsv")){
			summary.writer.println("panel\tseed\tmargin\tarm\tend\tmatched\texact\twithin3\twithin10\texact_pct\twithin3_pct\tboth_exact\tboth_within3\trms_bp\traw_rms_bp\trms_minus_raw_bp\tcomparator\tcomparator_rms_bp\trms_minus_comparator_bp\tbetter\tworse\ttied\tbetter_fraction\tcomparator_better\tcomparator_worse\tcomparator_tied\tmoved\tunmoved\texact_harmed\tnear_harmed\tdistant_improved\tdistant_worsened\tjoint_invalid\tcrossed\toutside_contig\tbelow_minimum\tabove_maximum\tcenter_held\tfixed\tmissing_features\tscored\tunsupported\tscope\tqualification");
			for(int seed=0; seed<labels.length; seed++){
				for(int scope=0; scope<2; scope++){
					for(int margin=0; margin<2; margin++){
						for(int arm=0; arm<4; arm++){
							for(int end=0; end<2; end++){
								if(stats[seed][scope][margin][arm][end]==null){continue;}
								b.clear().append(panel).tab().append(labels[seed]).tab().append(MARGINS[margin], 1).tab().append(ARMS[arm]).tab().append(ENDS[end]);
								stats[seed][scope][margin][arm][end].row(b, arm, scope);
								summary.writer.print(b);
							}
						}
					}
				}
			}
		}
	}

	public static int[] proposal(final int start, final int stop, final float[][] scores, final float margin){
		final int[] result=new int[2];
		proposal(start, stop, scores, margin, result);
		return result;
	}
	private static void proposal(final int start, final int stop, final float[][] scores, final float margin, final int[] result){
		require(scores!=null && scores.length==2 && result!=null && result.length==2, "Both ends share one unmoved candidate span");
		result[0]=Math.addExact(start, scores[0]==null ? 0 : NcrnaEndpointConfidence.selectFull(scores[0], 10, margin)-10);
		result[1]=Math.addExact(stop, scores[1]==null ? 0 : NcrnaEndpointConfidence.selectFull(scores[1], 10, margin)-10);
	}

	private boolean available(final int seed, final int end){return !netPaths[1][2*seed+end].equals("none");}

	private void writeTrial(final Output out, final ByteBuilder b, final Call call, final int seed, final int margin,
			final int arm, final int[] proposed, final int bad, final int[] ends, final int scope, final int target){
		b.clear().append(call.row).tab().append(call.input).tab().append(call.reference).tab().append(call.modelName).tab().append(call.seqid)
			.tab().append(call.minus ? "-" : "+").tab().append(call.length).tab().append(call.fixed).tab().append(support[call.model])
			.tab().append(call.start).tab().append(call.stop).tab().append(call.truthStart).tab().append(call.truthStop)
			.tab().append(call.proposedStart).tab().append(call.proposedStop).tab().append(GEOMETRY[call.centerBad])
			.tab().append(labels[seed]).tab().append(MARGINS[margin], 1).tab().append(ARMS[arm])
			.tab().append(proposed[0]).tab().append(proposed[1]).tab().append(GEOMETRY[bad])
			.tab().append(ends[0]).tab().append(ends[1]).tab().append((long)ends[0]-call.truthStart).tab().append((long)ends[1]-call.truthStop)
			.tab().append(SCOPES[scope]).tab().append(scope==0 ? ENDS[target] : "both").nl();
		out.writer.print(b);
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/

	private final class Stats{
		void add(final long error, final long raw, final long comparison, final long otherError,
				final int bad, final Call call, final int arm, final int end){
			n++;
			final long absolute=Math.abs(error);
			final long original=Math.abs(raw);
			if(absolute==0){exact++;}
			if(absolute<=3){within3++;}
			if(absolute<=10){within10++;}
			if(absolute==0 && otherError==0){bothExact++;}
			if(absolute<=3 && Math.abs(otherError)<=3){bothNear++;}
			if(absolute<original){better++;}else if(absolute>original){worse++;}else{tied++;}
			if(absolute<Math.abs(comparison)){compareBetter++;}else if(absolute>Math.abs(comparison)){compareWorse++;}else{compareTied++;}
			if(error!=raw){moved++;}
			if(original==0 && absolute>0){exactHarmed++;}
			if(original<=3 && absolute>3){nearHarmed++;}
			if(original>10 && absolute<original){distantImproved++;}
			if(original>10 && absolute>original){distantWorsened++;}
			squared+=(double)error*error;
			rawSquared+=(double)raw*raw;
			compareSquared+=(double)comparison*comparison;
			invalid[bad]++;
			if(call.centerBad!=0){centerHeld++;}
			if(support[call.model]<10){unsupported++;}
			if((call.fixed&(1<<end))!=0){fixed++;}
			else if(arm==1 || arm==3){
				final int view=arm==1 ? 0 : 1;
				if(call.features[view][end]!=null){scored++;}
				else if(view==0 || call.centerBad==0){missing++;}
			}
			assert(better+worse+tied==n && compareBetter+compareWorse+compareTied==n)
				: "Both paired comparisons must retain the complete frozen matched-call denominator";
		}
		/** Single-end rows omit pair statistics; no-change cohorts omit their movement fraction.
		 * Scope and unchanged population counts retain the corresponding presence conditions. */
		void row(final ByteBuilder b, final int arm, final int scope){
			require(n==expected, "Each arm must retain every matched real call");
			final double rms=Math.sqrt(squared/n);
			final double rawRms=Math.sqrt(rawSquared/n);
			final double compareRms=Math.sqrt(compareSquared/n);
			require(Double.isFinite(rms) && Double.isFinite(rawRms) && Double.isFinite(compareRms),
				"Computed real-call RMS values must be finite before report serialization");
			b.tab().append(n).tab().append(exact).tab().append(within3).tab().append(within10)
				.tab().append(100.0*exact/n, 8).tab().append(100.0*within3/n, 8);
			if(scope==0){b.append("\t\t");}else{b.tab().append(bothExact).tab().append(bothNear);}
			b.tab().append(rms, 8).tab().append(rawRms, 8).tab().append(rms-rawRms, 8).tab().append(arm==3 ? ARMS[1] : ARMS[0])
				.tab().append(compareRms, 8).tab().append(rms-compareRms, 8).tab().append(better).tab().append(worse).tab().append(tied).tab();
			if(better+worse>0){b.append((double)better/(better+worse), 8);}
			b.tab().append(compareBetter).tab().append(compareWorse).tab().append(compareTied).tab().append(moved).tab().append(n-moved)
				.tab().append(exactHarmed).tab().append(nearHarmed).tab().append(distantImproved).tab().append(distantWorsened)
				.tab().append(n-invalid[0]).tab().append(invalid[1]).tab().append(invalid[2]).tab().append(invalid[3]).tab().append(invalid[4])
				.tab().append(centerHeld).tab().append(fixed).tab().append(missing).tab().append(scored).tab().append(unsupported).tab().append(SCOPES[scope])
				.append("\tCONDITIONAL_MATCHED_REAL_CALLS_NOT_CALLER_REPRODUCTION").nl();
		}
		long n, exact, within3, within10, bothExact, bothNear, better, worse, tied, compareBetter, compareWorse, compareTied;
		long moved, exactHarmed, nearHarmed, distantImproved, distantWorsened, centerHeld, fixed, missing, scored, unsupported;
		double squared, rawSquared, compareSquared;
		final long[] invalid=new long[5];
	}

	private static final class Source{
		Source(final String id_, final String fasta_, final String assembly_){
			require(!assembly_.isEmpty() && !assembly_.equals(".") && !assembly_.equals("NA"),
				"Each real-call source must name its exact assembly for fold exclusion");
			id=id_;
			fasta=fasta_;
			assembly=assembly_;
		}
		final String id, fasta;
		String assembly;
		boolean canonicalAssembly;
		int ordinal=-1, savedOrdinal=-1;
		final HashMap<String, ArrayList<Call>> calls=new HashMap<String, ArrayList<Call>>();
		final HashMap<String, Call> keys=new HashMap<String, Call>();
	}
	private static final class Call{
		String key(){return seqid+'\t'+(minus ? "-" : "+")+'\t'+low+'\t'+high;}
		String status(final int arm, final int end){
			if((fixed&(1<<end))!=0){return "FIXED";}
			if(arm==1 && centerBad!=0){return "INVALID_CENTER_RAW_HOLD";}
			return features[arm][end]==null ? "MISSING_FEATURE_HOLD" : "AVAILABLE";
		}
		String input, assembly, reference, modelName, seqid;
		int model, ordinal, row, low, high, truthLow, truthHigh, length, expectedLength, start, stop, truthStart, truthStop;
		int fixed=-1, centerBad, centerStart, centerStop;
		long proposedStart, proposedStop;
		boolean minus;
		final int[] usable=new int[2];
		float[][][] features;
		/** No centered-span measurement until a gate feature end is usable; valid fractions are[0,1]. */
		float centerGC=-1;
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------       Fields / Constants     ----------------*/
	/*--------------------------------------------------------------*/

	private String in, out, sources, sourceFormat, fastaList, gff, familyRoot, ref, hbm, seeds, offsets, panel, assemblyAudit, foldMap;
	private String[] labels;
	private final String[][] tablePaths=new String[2][], netPaths=new String[2][];
	private int expected=-1, expectedWitnesses=-1, expectedInputs=-1;
	private int heldFold=-2;
	private final long[] foldCounts=new long[6];
	private NcrnaFamily family;
	private ObjectIntMap<String> modelIndex;
	private float[] plateaus;
	private final RrnaPositionalKmerTable.Table[][][] geometryTables=new RrnaPositionalKmerTable.Table[2][][];
	private int[] median5, median3;
	private long[] support;
	private final HashMap<String, Source> byInput=new HashMap<String, Source>();
	private final ArrayList<Source> orderedSources=new ArrayList<Source>();
	private final ArrayList<Call> calls=new ArrayList<Call>();
	private String familyName="cobalamin_rf00174", member;
	private boolean geometryOnly, gateVectors;
	public static final String LEGACY_GEOMETRY_HEADER="row1\tinput_id\tassembly\treference_id\temitted_member\tmodel\tseqid\tstrand\tcontig_length\tfixed_mask\traw_start0\traw_stop0\ttrue_start0\ttrue_stop0\tmin_length\tmax_length";
	/** The last cell is empty when the family has no maximum length. */
	public static final String GEOMETRY_HEADER=LEGACY_GEOMETRY_HEADER.replace("\tmax_length", "\tmax_length_if_bounded");
	public static final String GATE_HEADER=GEOMETRY_HEADER+"\tcenter_start0\tcenter_stop0\ttraining_n\tcenter_geometry\tend\tstatus\tgate_row1\tproposal0\tneural_delta_bp\tapplied_median_bp\tstrictly_better\tsingle_geometry";
	public static final String LEGACY_GATE_HEADER=GATE_HEADER.replace("\tmax_length_if_bounded", "\tmax_length");
	private static final String[] GATE_SUFFIXES={".vectors.tsv.gz", ".rows.tsv.gz", ".scores.tsv.gz", ".summary.tsv"};
	private static final String[] ENDS={"5prime", "3prime"};
	private static final String[] ARMS={"raw", "uncentered_nn", "center_only", "centered_nn"};
	private static final String[] SCOPES={"SINGLE_END_OTHER_AT_BASE", "PAIRED_JOINT"};
	private static final String[] GEOMETRY={"VALID", "CROSSED", "OUTSIDE_CONTIG", "BELOW_MINIMUM", "ABOVE_MAXIMUM"};
	private static final float[] MARGINS={0, 0.1f};
	private static final String[] SUFFIXES={".features.tsv", ".scores.tsv", ".trials.tsv", ".summary.tsv"};
	private static final String WITNESS_HEADER="input_id\treference_id\temitted_member\tmodel\tseqid\tcall_start\tcall_stop\tref_start\tref_stop\tstrand\td5\td3";
	private static final String SOURCE_HEADER="input_ordinal0\tassembly_id\tinput_id\tfasta\tcanonical_new_gff";
	private static final String LEGACY_HEADER="index\tgenome\tassembly\ttaxid\tdomain\tphylum\tgff\tfasta\tloci\tgenome_bases\tcontigs";
}
