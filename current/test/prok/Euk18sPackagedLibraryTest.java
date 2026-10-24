package test.prok;

import prok.*;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

/** Loads the actual default and explicit merged control through the production
 * entry point, including exact sequence-bound seed positions and cutoff rows.
 * @author Brian Bushnell, Raiden */
public final class Euk18sPackagedLibraryTest {
	public static void main(String[] args){
		if(args.length>1){throw new IllegalArgumentException("Optional historical seven-model control directory");}
		final Euk18sRuntimeConfig defaults=new Euk18sRuntimeConfig();defaults.parse("euk18s", "t");
		final NcrnaFamily nine=load(defaults);
		check(nine.library.length==9 && nine.joinedPositions.modelCount()==9, "Bare default loading must consume the selected nine models and their bound coordinates");
		for(int i=0; i<9; i++){
			check(!nine.modelNames[i].equals("euk18S_mitochondria") && nine.modelThresholds.pass(i)==.68f,
				"Removing the mitochondrial model must preserve every retained model's0.68 cutoff");
		}
		final Euk18sRuntimeConfig explicit=new Euk18sRuntimeConfig();explicit.parse("euk18s", "t");
		final FamilyResourceResolver resolver=new FamilyResourceResolver();
		explicit.parse("euk18sconsensus",resolver.resolve("euk18S","consensus",null).path);
		explicit.parse("euk18smodelcutoffs",resolver.resolve("euk18S","cutoffs",null).path);
		explicit.parse("euk18slivepositions",resolver.resolve("euk18S","joined_model_positions",null).path);
		final NcrnaFamily bound=load(explicit);
		check(Arrays.equals(nine.modelNames,bound.modelNames) && nine.joinedPositions.modelCount()==bound.joinedPositions.modelCount()
			&& nine.kmerSet.size()==bound.kmerSet.size() && nine.indexTopN==bound.indexTopN && nine.indexScoreMargin==bound.indexScoreMargin
			&& nine.windowPad==bound.windowPad && nine.pacBioCosts.sub==bound.pacBioCosts.sub,
			"Bare package and explicit selected resource bindings must resolve identical models and calling controls");
		for(int i=0;i<nine.library.length;i++){check(Arrays.equals(nine.library[i],bound.library[i]),"Default and explicit model sequences must be byte-identical");}
		if(args.length==0){System.out.println("EUK18S_PACKAGED_LIBRARY_PASS default9 explicit_default_parity");return;}
		final Euk18sRuntimeConfig control=new Euk18sRuntimeConfig();control.parse("euk18s", "t");
		control.parse("euk18sconsensus", new File(args[0], "consensus.fa").toString());
		control.parse("euk18smodelcutoffs", new File(args[0], "cutoffs.tsv").toString());
		control.parse("euk18slivepositions", new File(args[0], "seed_positions.tsv").toString());
		final NcrnaFamily seven=load(control);
		check(seven.library.length==7 && seven.joinedPositions.modelCount()==7, "Explicit original-v2 control must consume its own seven-model metadata");
		check(seven.seedMinHits==nine.seedMinHits && seven.seedDistinct==nine.seedDistinct && seven.indexK==nine.indexK
			&& seven.windowPad==nine.windowPad && seven.indexScoreMargin==nine.indexScoreMargin && seven.pacBioCosts.sub==nine.pacBioCosts.sub,
			"The library comparison must not silently change seed, shortlist, window or alignment settings");
		final Euk18sRuntimeConfig wrong=new Euk18sRuntimeConfig();wrong.parse("euk18s", "t");
		wrong.parse("euk18slivepositions", new File(args[0], "seed_positions.tsv").toString());
		boolean rejected=false;try{load(wrong);}catch(IllegalArgumentException expected){rejected=true;}
		check(rejected, "Default nine-model sequences must reject stale seven-model coordinates");
		System.out.println("EUK18S_PACKAGED_LIBRARY_PASS default9 explicit_control7 stale_positions_rejected unchanged_controls");
	}
	private static NcrnaFamily load(Euk18sRuntimeConfig config){
		final ArrayList<NcrnaFamily> list=new ArrayList<NcrnaFamily>();config.apply(list);
		check(list.size()==1, "Each isolated registration must produce exactly one complete family");return list.get(0);
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}
