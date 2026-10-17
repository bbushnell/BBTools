package prok;

import ml.CellNet;

/** Independent11/14-input tRNA adapter; loading has no process-wide side effects.
 * It never uses the generic10-input ncRNA feature vector or a mean-length feature.
 * @author Brian, Raiden */
final class TrnaBoundaryResources {
	TrnaBoundaryResources(String net5Path, String net3Path, String startPath, String stopPath, boolean shared){
		assert(net5Path!=null && net3Path!=null && startPath!=null && stopPath!=null)
			: "The resolver must supply complete tRNA boundary bindings before template publication";
		net5=TrnaBoundaryScorer.load(net5Path);
		net3=shared ? net5 : TrnaBoundaryScorer.load(net3Path);
		final int inputs=net5.numInputs();
		FamilyResourceManifest.require(inputs==11 || inputs==14,
			"tRNA boundary nets require11 or14 inputs from TrnaBoundaryVectorGen: "+net5Path+" has "+inputs);
		FamilyResourceManifest.require(net3.numInputs()==inputs,
			"Both tRNA boundary nets must use the same feature geometry: "+net5Path+" and "+net3Path);
		start=TrnaNinemerTableBuilder.loadTable(startPath);
		stop=TrnaNinemerTableBuilder.loadTable(stopPath);
		FamilyResourceManifest.require(start.type==TrnaBoundaryFeatures.BoundaryType.START,
			"tRNA start table has the wrong boundary type: "+startPath);
		FamilyResourceManifest.require(stop.type==TrnaBoundaryFeatures.BoundaryType.STOP,
			"tRNA stop table has the wrong boundary type: "+stopPath);
	}
	final CellNet net5, net3;
	final TrnaNinemerTableBuilder.LoadedTable start, stop;
}
