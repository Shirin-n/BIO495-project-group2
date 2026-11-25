import Foundation
import Vision
import AppKit
let root=URL(fileURLWithPath:CommandLine.arguments[1])
var output:[[String:Any]]=[]
let en=FileManager.default.enumerator(at:root, includingPropertiesForKeys:nil)!
for case let url as URL in en where url.lastPathComponent.hasSuffix("_busco.jpeg") {
 let req=VNRecognizeTextRequest(); req.recognitionLevel = .accurate; req.usesLanguageCorrection=false
 try VNImageRequestHandler(url:url).perform([req])
 let lines=(req.results ?? []).compactMap { $0.topCandidates(1).first?.string }
 output.append(["file":url.path,"lines":lines])
}
let d=try JSONSerialization.data(withJSONObject:output,options:[.prettyPrinted,.sortedKeys])
FileHandle.standardOutput.write(d)
