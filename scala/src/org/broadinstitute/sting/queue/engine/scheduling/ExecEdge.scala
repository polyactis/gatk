package org.broadinstitute.sting.queue.engine.scheduling

import org.apache.commons.lang.text.{StrLookup, StrSubstitutor}
import java.lang.String

class ExecEdge(private val templateCommandString: String)
        extends ResourceEdge {
  private var convertedCommandString: String = _
  def commandString = convertedCommandString

  override def traverse(graph: JobScheduler) = {
    // Lookup any variable using the target node, or any of it's input nodes.
    val sub = new StrSubstitutor(new StrLookup { def lookup(key: String) = graph.lookup(ExecEdge.this, key, null) })
    convertedCommandString = sub.replace(templateCommandString)
    graph.processExec(this)
  }
}