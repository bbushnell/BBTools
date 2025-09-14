package server;

import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;

import shared.Parse;
import shared.Timer;

public class StressTest {

	public static void main(String[] args) {
		long iterations=(args.length>0 ? Parse.parseKMG(args[0]) : 10000);
		String payload=(args.length>1 ? args[1] : "NC_012345");
		String message="https://taxonomy.jgi.doe.gov/accession/"+payload;
		Timer t=new Timer();
		HttpClient client = HttpClient.newBuilder()
				.version(HttpClient.Version.HTTP_2)
				.build();

		// Repeated queries to simulate your batching
		for(long i = 0; i < iterations; i++) {
			HttpRequest request = HttpRequest.newBuilder().uri(URI.create(message)).build();
			client.sendAsync(request, HttpResponse.BodyHandlers.ofString());
		}
		t.stop("Sent "+iterations+" queries in ");
	}
	
}
