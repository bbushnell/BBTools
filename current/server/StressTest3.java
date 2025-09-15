package server;

import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;

import shared.Parse;
import shared.Timer;

public class StressTest3 {
    public static void main(String[] args) {
        System.setProperty("javax.net.debug", "ssl:handshake");
        
        long iterations = (args.length > 0 ? Parse.parseKMG(args[0]) : 10000);
        String payload = (args.length > 1 ? args[1] : "NC_012345");
        String message = "https://taxonomy.jgi.doe.gov/accession/" + payload;
        
        HttpClient client = HttpClient.newBuilder()
            .version(HttpClient.Version.HTTP_2)
            .build();
        
        Timer t = new Timer();
        for (long i = 0; i < iterations; i++) {
        	HttpResponse<String> response=null;
            try {
                HttpRequest request = HttpRequest.newBuilder().uri(URI.create(message)).build();
                //System.out.println("Sending request " + i);
                response = client.send(request, HttpResponse.BodyHandlers.ofString());
                //System.out.println("Response status: " + response.statusCode());
            } catch (Exception e) {
                System.err.println("Error on request " + i + ": " + e);
                e.printStackTrace();
            }
            assert(response!=null);
        }
        
        t.stop("Sent " + iterations + " queries in ");
    }
}
