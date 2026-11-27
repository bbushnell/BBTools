package prok;

import java.io.File;
import java.io.IOException;

import dna.AminoAcid;
import dna.GeneticCode;
import parse.LineParser1;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Tools;
import structures.ByteBuilder;

/** Prepared capture identity plus the actual loaded fixed-model runtime.
 * The offline worker verifies immutable resource provenance once per distinct
 * payload. This class neither hashes files nor substitutes names for that check.
 * A later JVM compares its freshly loaded native metadata to this exact record.
 * @author Brian, Raiden
 */
public final class NcrnaPathBindings{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public NcrnaPathBindings(final String provenance, final String recipe,
			final int sources, final GeneCaller caller, final String[] modelSpecs,
			final String[] profiles){
		require(caller!=null && modelSpecs!=null && modelSpecs.length>0 && sources>0,
			"Bindings require the actual loaded caller, ordered PGM specifications and declared source population");
		require("shipped".equals(recipe) || "working".equals(recipe), "Bind one explicit comparison recipe");
		require(recipe.equals("working") ? profiles!=null && profiles.length==3 : profiles==null,
			"The working recipe binds all three external ordering profiles; shipped binds none");
		caller.verifyPathReplayPolicy();
		merged=recipe.equals("working");
		data.append("#NCRNA_PATH_BINDINGS\t1\n");
		provenance(provenance);
		text("RECIPE", recipe);
		data.append("SOURCES\t").append(sources).nl();
		data.append("PGM_MERGE\t").append(PGMTools.normalize ? 1 : 0).tab()
			.append(GeneModel.ADD_CDS_ONLY ? 1 : 0).nl();
		for(int i=0; i<modelSpecs.length; i++){
			final String spec=modelSpecs[i];
			require(spec!=null && !spec.isEmpty(), "Every ordered weighted PGM specification must be explicit");
			data.append("PGM\t").append(i).tab();
			NcrnaPathGraphIO.text(data, spec).tab();
			final int weight=spec.indexOf('@');
			NcrnaPathGraphIO.text(data, realPath(weight<0 ? spec : spec.substring(0, weight))).nl();
		}
		path("COBALAMIN_ROOT", CallGenes.COBALAMIN_RESOURCE_ROOT);
		path("FAMILY_ROOT", CallGenes.pathFamilyRoot());
		path("CONSENSUS_OVERRIDE", CallGenes.NCRNA_CONSENSUS_OVERRIDE);
		path("MODELS_OVERRIDE", CallGenes.NCRNA_MODELS_OVERRIDE);
		path("SEEDS_OVERRIDE", CallGenes.NCRNA_KMERS_OVERRIDE);
		if(profiles!=null){
			for(int i=0; i<profiles.length; i++){
				require(profiles[i]!=null && !profiles[i].isEmpty(), "Each bound ordering profile needs an actual resource path");
				path("PROFILE_"+i, profiles[i]);
			}
		}
		geometry(caller.pgm);
		code(caller.pathGeneticCode());
		families();
		caller.pathPolicy().append(data).nl();
		data.append("COMPLETE\n");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** The worker owns the literal identities and verifies payloads before loading. */
	private void provenance(final String path){
		final LineParser1 p=new LineParser1('\t');
		try(Input input=new Input(path)){
			RrnaResourceIO.header(input.reader, "#NCRNA_PATH_PROVENANCE\t1", path);
			for(int i=0; i<PROVENANCE.length; i++){
				final byte[] row=input.reader.nextLine();
				require(row!=null, "Capture provenance is missing a required source/resource/configuration identity");
				p.set(row);
				require(p.terms()==2 && p.termEquals(PROVENANCE[i], 0), "Capture provenance fields must use their declared order and schema");
				final String value=NcrnaPathGraphIO.string(p, 1);
				if(i>0){require(hex(value, i==4 ? 40 : 20), "Provenance requires sha80 identities, except the Git commit identifier");}
				text(PROVENANCE[i], value);
			}
			require(input.reader.nextLine()==null, "Capture provenance cannot contain undeclared identities");
		}
	}
	private void geometry(final GeneModel model){
		assert(model!=null) : "The binding describes the actual loaded model used by native overlap rescoring";
		for(int i=0; i<model.allContainers.length; i++){
			final StatsContainer container=model.allContainers[i];
			final FrameStats[] parts=container.statsArray;
			for(int j=0; j<parts.length; j++){
				final FrameStats part=parts[j];
				require(part!=null, "Effective model geometry must be initialized before binding");
				data.append("GEOMETRY\t").append(i).tab().append(j).tab().append(part.k)
					.tab().append(part.frames).tab().append(part.leftOffset).nl();
			}
		}
	}
	private void code(final GeneticCode code){
		data.append("CODE\t").append(code==null ? "legacy" : "explicit").tab().append(code==null ? -1 : code.id()).nl();
		for(int i=0; i<64; i++){
			final int amino=code==null ? AminoAcid.codeToByte[i] : code.aminoAcid(i);
			final boolean start=code==null ? GeneModel.isStartCodon(i) : code.isStart(i);
			final boolean stop=code==null ? GeneModel.isStopCodon(i) : code.isStop(i);
			data.append("CODON\t").append(i).tab().append(amino).tab().append(start ? 1 : 0)
				.tab().append(stop ? 1 : 0).nl();
		}
	}
	private void families(){
		data.append("FAMILIES\t").append(GeneCaller.ncrnaFamilies.size()).nl();
		for(int i=0; i<GeneCaller.ncrnaFamilies.size(); i++){
			final NcrnaFamily family=GeneCaller.ncrnaFamilies.get(i);
			for(int member=0; member<MEMBERS.length; member++){
				if(MEMBERS[member].equals(family.name)){
					require(members[member]==null, "The bound clan axis contains one resource per registered member");
					members[member]=family;
				}
			}
			require(family.modelNames!=null && family.library!=null && family.modelNames.length==family.library.length,
				"The binding needs the actual registered index axis, not an inferred model-name list");
			data.append("FAMILY\t").append(i).tab();
			NcrnaPathGraphIO.text(data, family.name).tab().append(family.modelNames.length).tab()
				.append(Float.floatToRawIntBits(family.scoreA)).tab().append(Float.floatToRawIntBits(family.scoreB)).nl();
			for(int j=0; j<family.modelNames.length; j++){
				data.append("MODEL\t").append(i).tab().append(j).tab();
				NcrnaPathGraphIO.text(data, family.modelNames[j]).tab().append(family.library[j].length).nl();
			}
		}
		for(final NcrnaFamily member:members){require(member!=null, "Bind all three actual member index axes before saved scoring");}
	}
	/** Native combined ordinals concatenate the three ordered member axes
	 * (NcrnaClanScavenger.Models); separate callers retain member-local ordinals. */
	public void verifyScoreBasis(final String family, final String model, final int nativeModel){
		require(nativeModel>=0, "A saved score basis needs its actual native model ordinal");
		int local=nativeModel;
		for(final NcrnaFamily member:members){
			if(merged){
				if(local<member.modelNames.length){
					require(member.name.equals(family) && member.modelNames[local].equals(model),
						"Saved score owner/model does not match the bound combined native index");
					return;
				}
				local-=member.modelNames.length;
			}else if(member.name.equals(family)){
				require(local<member.modelNames.length && member.modelNames[local].equals(model),
					"Saved score owner/model does not match the bound member-local native index");
				return;
			}
		}
		throw new IllegalArgumentException("Saved score basis lies outside the bound target model axes");
	}
	private void text(final String key, final String value){
		data.append(key).tab();
		NcrnaPathGraphIO.text(data, value).nl();
	}
	private void path(final String key, final String value){
		data.append(key).tab();
		NcrnaPathGraphIO.text(data, value).tab();
		if(value!=null){NcrnaPathGraphIO.text(data, realPath(value));}
		data.nl();
	}
	public void write(final String path){
		require(Tools.testOutputFiles(false, false, false, path), "Prepared binding output must be fresh");
		try(Output output=new Output(path)){output.writer.print(data);}
	}
	/** Compare semantic native rows without a second resource hash or genome read. */
	public void verify(final String path){
		int cursor=0, line=0;
		try(Input input=new Input(path)){
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				line++;
				for(final byte value:row){
					require(cursor<data.length() && data.array[cursor++]==value,
						"Fresh native metadata differs from the prepared capture binding at line "+line);
				}
				require(cursor<data.length() && data.array[cursor++]=='\n', "Prepared binding row is truncated or has extra fields at line "+line);
			}
		}
		require(cursor==data.length(), "Prepared binding ended before the complete native runtime description");
	}
	/** Graph framing names the full prepared binding, not an independently replaceable manifest. */
	public static String graphHeader(final String manifestSha80){
		require(hex(manifestSha80, 20), "Bound graph headers require the prepared manifest's sha80");
		return "#NCRNA_PATH_GRAPHS\t2\t"+manifestSha80;
	}
	private static String realPath(final String path){
		try{return new File(path).getCanonicalPath();}
		catch(IOException failure){throw new IllegalArgumentException("Cannot resolve a bound resource path: "+path, failure);}
	}
	private static boolean hex(final String value, final int length){
		if(value==null || value.length()!=length){return false;}
		for(int i=0; i<value.length(); i++){
			final char c=value.charAt(i);
			if(!(c>='0' && c<='9' || c>='a' && c<='f')){return false;}
		}
		return true;
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}
	private final ByteBuilder data=new ByteBuilder();
	private final boolean merged;
	private final NcrnaFamily[] members=new NcrnaFamily[3];
	private static final String[] MEMBERS={"cobalamin_rf00174", "cobalamin_rf01482", "cobalamin_rf01689"};
	private static final String[] PROVENANCE={"capture", "source_roster_sha80", "resource_catalog_sha80",
		"caller_config_sha80", "code_commit", "class_seal"};
}
