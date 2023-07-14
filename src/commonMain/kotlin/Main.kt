import korlibs.korge.*
import korlibs.korge.scene.*
import korlibs.korge.view.*
import org.korge.samples.mymodule.*

suspend fun main() = Korge {
    sceneContainer().changeTo({ MainMyModuleScene() })
}

class MainMyModuleScene : Scene() {
    override suspend fun SContainer.sceneMain() {
        text(MyModule.TEXT)
    }
}